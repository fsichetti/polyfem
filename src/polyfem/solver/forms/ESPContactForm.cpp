#include "ESPContactForm.hpp"
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/Types.hpp>
#include <polyfem/utils/Timer.hpp>
#include <polyfem/utils/MatrixUtils.hpp>
#include <polyfem/utils/MaybeParallelFor.hpp>
#include <polyfem/quadrature/TriQuadrature.hpp>

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>
#include <ipc/barrier/barrier.hpp>

namespace polyfem::solver
{
	namespace
	{
		ipc::FaceQuadRule build_quad_rule(const int order)
		{
			polyfem::quadrature::Quadrature quad;
			polyfem::quadrature::TriQuadrature(true).get_quadrature(order, quad);

			ipc::FaceQuadRule rule;
			rule.reserve(quad.size());
			for (int i = 0; i < quad.size(); i++)
			{
				const double x = quad.points(i, 0);
				const double y = quad.points(i, 1);
				// polyfem divides weights by 2 (reference triangle area); multiply
				// back to match ipc-toolkit's convention of weights summing to 1.
				const double w = quad.weights(i) * 2.0;
				rule.push_back({{{1.0 - x - y, x, y}}, w});
			}
			return rule;
		}
	} // namespace

	/// Check if a quadrature point is at a triangle vertex (one barycentric
	/// coordinate ≈ 1, others ≈ 0).
	bool is_vertex_point(const ipc::FaceQuadPoint &qp, int vertex, double tol = 1e-12)
	{
		for (int i = 0; i < 3; i++)
		{
			if (std::abs(qp.lambda[i] - (i == vertex ? 1.0 : 0.0)) > tol)
				return false;
		}
		return true;
	}

	/// Verify the quadrature rule contains the 3 triangle vertices.
	void verify_vertices_in_quad_rule(const ipc::FaceQuadRule &rule)
	{
		for (int v = 0; v < 3; v++)
		{
			bool found = false;
			for (const auto &qp : rule)
			{
				if (is_vertex_point(qp, v))
				{
					found = true;
					break;
				}
			}
			if (!found)
				throw std::runtime_error("Face quadrature rule is missing triangle vertex " + std::to_string(v) + "; choose a quadrature scheme that includes corner points");
		}
	}

	ipc::EEMollifierType ee_mollifier_type(const std::string &name)
	{
		if (name == "smoothstep_rel")
			return ipc::EEMollifierType::SMOOTHSTEP_REL;
		if (name == "smoothstep")
			return ipc::EEMollifierType::SMOOTHSTEP;
		if (name == "smoothstep_sqrt")
			return ipc::EEMollifierType::SMOOTHSTEP_SQRT;
		if (name == "gcp")
			return ipc::EEMollifierType::GCP;
		if (name == "gcp_fixed")
			return ipc::EEMollifierType::GCP_FIXED;
		if (name == "gcp_sqrt")
			return ipc::EEMollifierType::GCP_SQRT;
		log_and_throw_error("Unknown ESP edge-edge mollifier \"{}\".", name);
	}

	ipc::ESPParameters init_params(const double dhat, const json &esp_params, const bool skip_obstacles, std::shared_ptr<ipc::Barrier> barrier, const int dim)
	{
		// quadrature_order is "auto" or an order >= 1.
		const json order = esp_params.value("quadrature_order", json("auto"));
		const bool auto_rule = order.is_string();
		if (auto_rule && order.get<std::string>() != "auto")
			log_and_throw_error("ESP quadrature_order must be \"auto\" or an integer >= 1, not \"{}\".", order.get<std::string>());
		const int quadrature_order = auto_rule ? 0 : order.get<int>();
		if (!auto_rule && quadrature_order < 1)
			log_and_throw_error("ESP quadrature_order must be \"auto\" or an integer >= 1, not {}.", quadrature_order);
		const double dbar_factor = esp_params["dbar_factor"];
		const bool area_weights = esp_params["area_weights"];
		const ipc::ESPParameters::IntegrationType itype = skip_obstacles ? ipc::ESPParameters::IntegrationType::NO_OBST : ipc::ESPParameters::IntegrationType::NORMAL;
		// In 2D the order is the Gauss-Lobatto edge order; "auto" means 1
		// there. 3D does not use it.
		const int edge_order = auto_rule ? 1 : quadrature_order;
		ipc::ESPParameters params(dhat, dbar_factor, edge_order, area_weights, itype);
		params.barrier = barrier ? barrier : (dim == 3 ? std::make_shared<ipc::InversePowerBarrier<>>(2.0) : std::make_shared<ipc::InversePowerBarrier<>>(1.0));
		if (!auto_rule)
		{
			params.set_quad_rule(build_quad_rule(quadrature_order));
			verify_vertices_in_quad_rule(params.quad_rule());
		}
		else if (dim == 3)
		{
			// Per-face composite vertex (lattice) rules chosen from the rest
			// face area, so the points are about lattice_spacing * dhat apart.
			params.set_auto_lattice(
				esp_params.value("lattice_spacing", 1.0) * dhat,
				esp_params.value("max_lattice_subdivisions", 16));
		}

		// Read with the input spec's defaults, for callers that pass a
		// partial object.
		params.legacy_element_weights = esp_params.value("legacy_element_weights", false);
		params.ee_mollifier.type = ee_mollifier_type(esp_params.value("ee_mollifier", "smoothstep_rel"));
		params.ee_mollifier.rel_eps = esp_params.value("ee_mollifier_eps", 1.0);
		params.ee_mollifier.width = esp_params.value("ee_mollifier_width", 150e-6);
		params.ee_mollifier.length = esp_params.value("ee_mollifier_length", 1.5e-3);
		params.ee_mollifier_power = esp_params.value("ee_mollifier_power", 0);
		params.ee_weight = esp_params.value("ee_weight", 1.5);
		params.ee_area_scaling = esp_params.value("ee_area_scaling", true);
		params.ee_area_ref = esp_params.value("ee_area_ref", 0.0);

		// Local transfer needs normalized weights and a face rule whose
		// corners carry weight (they donate it).
		bool transfer = esp_params.value("local_transfer", true);
		if (transfer && !esp_params.value("normalize_weights", true))
		{
			logger().info("ESP local_transfer needs normalize_weights; turning it off.");
			transfer = false;
		}
		const auto &cw = params.corner_weights();
		if (transfer && params.lattice_spacing() <= 0 && !(cw[0] > 0 && cw[1] > 0 && cw[2] > 0))
		{
			logger().info("ESP local_transfer needs a face rule with positive corner weights; turning it off.");
			transfer = false;
		}
		params.local_quadrature_transfer = transfer;
		return params;
	}

	ESPContactForm::ESPContactForm(const ipc::CollisionMesh &collision_mesh,
								   const double dhat,
								   const double avg_mass,
								   const json esp_params,
								   const bool skip_obstacles,
								   std::shared_ptr<ipc::Barrier> barrier,
								   const bool use_adaptive_dhat,
								   const bool use_adaptive_barrier_stiffness,
								   const bool is_time_dependent,
								   const bool enable_shape_derivatives,
								   const ipc::BroadPhaseMethod broad_phase_method,
								   const double ccd_tolerance,
								   const int ccd_max_iterations,
								   const double dhat_epsilon_scale) : ContactForm(collision_mesh, dhat, avg_mass, use_adaptive_barrier_stiffness, is_time_dependent, enable_shape_derivatives, broad_phase_method, ccd_tolerance, ccd_max_iterations, dhat_epsilon_scale), params(init_params(dhat, esp_params, skip_obstacles, barrier, static_cast<int>(collision_mesh.dim()))),
																	  barrier_potential_(params, esp_params["normalize_weights"])
	{
		if (use_adaptive_dhat)
			log_and_throw_error("Adaptive dhat is not implemented for ESPContactForm!");
	}

	void ESPContactForm::update_barrier_stiffness(const Eigen::VectorXd &x, const Eigen::MatrixXd &grad_energy)
	{
		if (!use_adaptive_barrier_stiffness())
			return;

		log_and_throw_error("Adaptive barrier stiffness not implemented for ESPContactForm!");
	}

	void ESPContactForm::update_collision_set(const Eigen::MatrixXd &displaced_surface)
	{
		// Store the previous value used to compute the constraint set to avoid duplicate computation.
		if (cached_displaced_surface.size() == displaced_surface.size() && cached_displaced_surface == displaced_surface)
			return;

		collision_set_.build(
			collision_mesh_, displaced_surface, params, broad_phase_.get());
		cached_displaced_surface = displaced_surface;
	}

	double ESPContactForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		const Eigen::MatrixXd displaced = compute_displaced_surface(x);
		if (cached_displaced_surface != displaced)
		{
			return 0.;
		}
		return barrier_potential_(collision_set_, collision_mesh_, displaced);
	}

	Eigen::VectorXd ESPContactForm::value_per_element_unweighted(const Eigen::VectorXd &x) const
	{
		log_and_throw_error("value_per_element_unweighted not implemented!");
	}

	void ESPContactForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		const Eigen::MatrixXd displaced = compute_displaced_surface(x);
		if (cached_displaced_surface != displaced)
		{
			gradv.setZero(x.size());
			return;
		}
		gradv = barrier_potential_.gradient(collision_set_, collision_mesh_, displaced);
		gradv = collision_mesh_.to_full_dof(gradv);
	}

	void ESPContactForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("barrier hessian");
		hessian = barrier_potential_.hessian(collision_set_, collision_mesh_, compute_displaced_surface(x), project_to_psd_ ? ipc::PSDProjectionMethod::CLAMP : ipc::PSDProjectionMethod::NONE);
		hessian = collision_mesh_.to_full_dof(hessian);
	}

	void ESPContactForm::post_step(const polysolve::nonlinear::PostStepData &data)
	{
		const Eigen::MatrixXd displaced_surface = compute_displaced_surface(data.x);

		// Always requires update_collision_set
		update_collision_set(displaced_surface);

		const double curr_distance = collision_set_.compute_minimum_distance(collision_mesh_, displaced_surface);
		if (!std::isinf(curr_distance))
		{
			const double ratio = sqrt(curr_distance) / dhat();
			const auto log_level = (ratio < 1e-6) ? spdlog::level::err : ((ratio < 1e-4) ? spdlog::level::warn : spdlog::level::debug);
			polyfem::logger().log(log_level, "Minimum distance during solve: {}, dhat: {}", sqrt(curr_distance), dhat());
		}

		if (data.iter_num == 0)
			return;

		if (use_adaptive_barrier_stiffness_)
		{
			if (is_time_dependent_)
			{
				const double prev_barrier_stiffness = barrier_stiffness();

				barrier_stiffness_ = ipc::update_barrier_stiffness(
					prev_distance_, curr_distance, max_barrier_stiffness_,
					barrier_stiffness(), ipc::world_bbox_diagonal_length(displaced_surface),
					dhat_epsilon_scale_);

				if (barrier_stiffness() != prev_barrier_stiffness)
				{
					polyfem::logger().debug(
						"updated barrier stiffness from {:g} to {:g}",
						prev_barrier_stiffness, barrier_stiffness());
				}
			}
			else
			{
				// TODO: missing feature
				// update_barrier_stiffness(data.x);
			}
		}

		prev_distance_ = curr_distance;
	}
} // namespace polyfem::solver