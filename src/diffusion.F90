module diffusion
    use mod_tool, only: cal_cfl_time_step
    use mod_vdiff, only: vdiff_by_k_theory
    use mod_hdiff, only: cal_kh_by_deformation_method, hdiff1d_by_k_theory
end module diffusion