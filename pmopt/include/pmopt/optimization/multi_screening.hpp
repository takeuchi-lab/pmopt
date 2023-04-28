#pragma once

namespace PMOpt
{
    namespace Optimization
    {

        /**
         * @brief compute squared distance of dual variables
         * @tparam Loss
         * @param vars_1
         * @param vars_2
         * @param param
         * @return double 
         */
        template < class Loss >
        double distance(
            const Variables & vars_1,
            const Variables & vars_2,
            const Parameters & params
        )
        noexcept
        {
            double distance = 0.0;
            for (auto i_instance : * params._instance_indices.get())
            {
                auto r_1 = Loss::dual(vars_1, params, i_instance) * vars_1._scale;
                auto r_2 = Loss::dual(vars_2, params, i_instance) * vars_2._scale;
                auto m = params._mask[i_instance];
                distance += (r_1 - r_2) * (r_1 - r_2) * m;
            }
            return std::sqrt(distance);
        }


        /**
         * @brief true if sphere of vars_1 is including sphere of vars_2
         * @param vars_1 
         * @param vars_2 
         * @param distance 
         * @return bool 
         */
        inline bool is_including(
            const Variables & vars_1,
            const Variables & vars_2,
            double distance
        )
        noexcept
        {
            // If '<' -> '<=', no sphere will be used
            // in the case that vars_1._radius == vars_2.radius.
            return distance < vars_1._radius - vars_2._radius;
        }


        /**
         * @brief check if optimal solution of vars_1 is in sphere of vars_2 
         * @param radius_1 
         * @param radius_2 
         * @param d_score 
         * @param norm 
         * @param distance 
         * @return bool 
         */
        inline bool is_in(
            double radius_1,
            double radius_2,
            double d_score,
            double norm,
            double distance
        )
        noexcept
        {
            auto d_r = radius_1 * radius_1 - radius_2 * radius_2;
            return distance * distance + 2 * radius_1 * d_score / norm + d_r <= 0;
        }


        /**
         * @brief check inclusion relationship between spheres
         * @note NOT compute distance_matrix[i][j] for all i, j
         * @tparam Loss 
         * @param[out] distance_matrix 
         * @param[out] including 
         * @param vars 
         * @param params
         * @param size 
         */
        template < class Loss >
        void check_inclusion(
            std::vector<std::vector<double>> & distance_matrix, 
            std::vector<unsigned int> & including,
            const std::vector<Variables> & vars,
            const Parameters & params,
            unsigned int n_vars
        )
        noexcept
        {
            // initialize vectors
            distance_matrix.resize(n_vars, std::vector<double>(n_vars));
            including.resize(n_vars, false);

            // compute distance matrix and inclusion relationship
            for (unsigned int i_vars_1 = 0; i_vars_1 < n_vars; ++i_vars_1)
            {
                // if this sphere includes others then it is inactive
                if (including[i_vars_1])
                    continue;

                for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < n_vars; ++i_vars_2)
                {
                    // if this sphere includes others then it is inactive
                    if (including[i_vars_2])
                        continue;

                    const auto & vars_1 = vars[i_vars_1];
                    const auto & vars_2 = vars[i_vars_2];
                    auto distance = Optimization::distance<Loss>(vars_1, vars_2, params);
                    distance_matrix[i_vars_1][i_vars_2] = distance;
                    including[i_vars_1] = is_including(vars_1, vars_2, distance);
                    including[i_vars_2] = is_including(vars_2, vars_1, distance);
                }
            }
        }
        


        /**
         * @brief compute screening and pruning score
         * @param[out] pruning_criterion
         * @param[out] screening_criterion 
         * @param distance 
         * @param scale_1
         * @param scale_2 
         * @param radius_1
         * @param radius_2 
         * @param dot_upper_1
         * @param dot_upper_2 
         * @param dot_lower_1 
         * @param dot_lower_2 
         * @param x_sum
         * @param x_sum_sq
         * @param n
         */
        void safe_pruning_multi(
            double & pruning_criterion,
            double & screening_criterion,
            double distance,
            double scale_1,
            double scale_2,
            double radius_1,
            double radius_2,
            double dot_upper_1,
            double dot_upper_2,
            double dot_lower_1,
            double dot_lower_2,
            double x_sum,
            double x_sum_sq,
            unsigned int n
        )
        noexcept
        {
            auto r_1 = radius_1;
            auto r_2 = radius_2;
            auto d = distance;

            // internal division ratio for center of intersecion sphere
            auto t = 0.5 + (r_2 * r_2 - r_1 * r_1) / d / d * 0.5;

            // radius of intersection sphere
            auto r_intersection = std::sqrt(r_2 * r_2 - t * t * d * d);

            // screening and pruning score of each sphere (already scaled)
            dot_upper_1 *= scale_1;
            dot_lower_1 *= scale_1;
            dot_upper_2 *= scale_2;
            dot_lower_2 *= scale_2;
            auto dot_value_1 = dot_upper_1 + dot_lower_1;
            auto dot_value_2 = dot_upper_2 + dot_lower_2;

            // difference of two screening scores
            auto d_dot = dot_value_1 - dot_value_2;

            // weighted mean of two dot products (positive / negative part respectively)
            auto dot_upper_intersection = dot_upper_1 * t + dot_upper_2 * (1 - t);
            auto dot_lower_intersection = dot_lower_1 * t + dot_lower_2 * (1 - t);
            auto dot_value_intersection = dot_upper_intersection + dot_lower_intersection;

            // TODO if feature is binary, tighter bound can be used
            // compute radial component of screening spheres
            auto radial_pruning = std::sqrt(x_sum_sq);
            auto radial_screening = std::sqrt(x_sum_sq - x_sum * x_sum / n);
            auto radial_intersection = std::sqrt(radial_screening * radial_screening - d_dot * d_dot / d / d);

            // screening criterion
            // check maximum value
            double screening_max;
            if (is_in(r_1, r_2, d_dot, radial_screening, d)) // if max of vars_1 is global max
                screening_max = dot_value_1 + radial_screening * r_1;
            else if (is_in(r_2, r_1, - d_dot, radial_screening, d)) // if max of vars_2 is global max
                screening_max = dot_value_2 + radial_screening * r_2;
            else // if max of inersection of vars_1 and vars_2 is global max
                screening_max = dot_value_intersection + radial_intersection * r_intersection;

            // check minimum value
            double screening_min;
            if (is_in(r_1, r_2, - d_dot, radial_screening, d)) // if min of vars_1 is global min
                screening_min = dot_value_1 - radial_screening * r_1;
            else if (is_in(r_2, r_1, d_dot, radial_screening, d)) // if min of vars_2 is global min
                screening_min = dot_value_2 - radial_screening * r_2;
            else // if min of intersection of vars_1 and vars_2 is global min
                screening_min = dot_value_intersection - radial_intersection * r_intersection;

            // max of abs is max of min and max
            // max |v| = max{max v, - min v}
            screening_criterion = std::max(screening_max, - screening_min);

            // pruning criterion
            // can not be calculated in the same way as screening
            auto dot_bound_1 = std::max(dot_upper_1, - dot_lower_1);
            auto dot_bound_2 = std::max(dot_upper_2, - dot_lower_1);
            pruning_criterion = std::min(dot_bound_1 + radial_pruning * r_1, dot_bound_2 + radial_pruning * r_2);
        }


        /**
         * @brief compute screening score woth multi reference
         * @param[out] criterion 
         * @param distance 
         * @param scale_1
         * @param scale_2 
         * @param radius_1
         * @param radius_2 
         * @param dot_value_1
         * @param dot_value_2 
         * @param x_sum
         * @param x_sum_sq
         * @param n
         */
        void safe_screening_multi(
            double & criterion,
            double distance,
            double scale_1,
            double scale_2,
            double radius_1,
            double radius_2,
            double dot_value_1,
            double dot_value_2,
            double x_sum,
            double x_sum_sq,
            unsigned int n
        )
        noexcept
        {
            auto r_1 = radius_1;
            auto r_2 = radius_2;
            auto d = distance;

            // internal division ratio for center of intersecion sphere
            auto t = 0.5 + (r_2 * r_2 - r_1 * r_1) / d / d * 0.5;

            // radius of intersection sphere
            auto r_intersection = std::sqrt(r_2 * r_2 - t * t * d * d);

            // difference of two screening scores
            dot_value_1 *= scale_1;
            dot_value_2 *= scale_2;
            auto d_dot = dot_value_1 - dot_value_2;

            // weight screening score
            auto dot_value_intersection = dot_value_1 * t + dot_value_2 * (1 - t);

            // compute radial component of screening spheres
            auto radial_screening = std::sqrt(x_sum_sq - x_sum * x_sum / n);
            auto radial_intersection = std::sqrt(radial_screening * radial_screening - d_dot * d_dot / d / d);

            // check maximum value
            double screening_max;
            if (is_in(r_1, r_2, d_dot, radial_screening, d)) // if max of vars_1 is global max
                screening_max = dot_value_1 + radial_screening * r_1;
            else if (is_in(r_2, r_1, - d_dot, radial_screening, d)) // if max of vars_2 is global max
                screening_max = dot_value_2 + radial_screening * r_2;
            else // if max of inersection of vars_1 and vars_2 is global max
                screening_max = dot_value_intersection + radial_intersection * r_intersection;

            // check minimum value
            double screening_min;
            if (is_in(r_1, r_2, - d_dot, radial_screening, d)) // if min of vars_1 is global min
                screening_min = dot_value_1 - radial_screening * r_1;
            else if (is_in(r_2, r_1, d_dot, radial_screening, d)) // if min of vars_2 is global min
                screening_min = dot_value_2 - radial_screening * r_2;
            else // if min of intersection of vars_1 and vars_2 is global min
                screening_min = dot_value_intersection - radial_intersection * r_intersection;

            // max of abs is max of min and max
            // max |v| = max{max v, - min v}
            criterion = std::max(screening_max, - screening_min);
        }
    }
}