function [q_s, q_a, q_p, q_sat0] = radiation(sigma, temp_s, temp_p, r_s, r_sp, r_sat, r_p, abs_sat, ref_p, emi_sat, vf_s_p, vf_sat_surr, vf_p_sat, emi_p )

q_s = (sigma * power(temp_s,4) * 4 * pi * power(r_s,2) * power(r_sat,2) * abs_sat) / (4 * power(r_sp,2));

q_a = (sigma * power(temp_s,4) * 4 * pi * power(r_s,2) * vf_s_p * ref_p * vf_p_sat * abs_sat);

q_p = (sigma * emi_p * power(temp_p,4) * 4 * pi * power(r_p,2) * vf_p_sat * abs_sat );

q_sat0 = (sigma * emi_sat * 4 * pi * power(r_sat,2) * vf_sat_surr );

end