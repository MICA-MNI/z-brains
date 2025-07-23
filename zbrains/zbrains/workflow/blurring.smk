sfwm_depths = config['sfwm_depths']

# --- Laplace Solver Rule ---
rule laplace_solver:
    input:
        aparc_aseg = "{input_dir}/{subject}_{session}_aparc+aseg.nii.gz"
    output:
        laplace = "{output_dir}/{subject}_{session}_laplace-wm.nii.gz"
    shell:
        """
        python -m sWM.laplace_solver {input.aparc_aseg} {output.laplace}
        """

# Use a wildcard {depth} for sfwm surfaces, expanded from config['sfwm_depths']

# --- Shift Surface Rule ---
rule shift_surface:
    input:
        white = "{input_dir}/{subject}_{session}_hemi-{hemi}_label-white.surf.gii",
        laplace = "{output_dir}/{subject}_{session}_laplace-wm.nii.gz"
    output:
        sfwm = "{output_dir}/{subject}/{session}/structural/{subject}_{session}_{hemi}_sfwm-{depth}mm.surf.gii"
    params:
        depths = lambda w: config['sfwm_depths']
    shell:
        """
        python -m sWM.surface_generator {input.white} {input.laplace} {output.sfwm} {wildcards.depth}
        """

# --- Volume to Surface Mapping: midthickness ---
rule v2s_midthickness:
    input:
        volumemap = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/maps/{{subject}}_{{session}}_space-nativepro_map-{{feature}}.nii.gz",
        surf = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/surf/{{subject}}_{{session}}_hemi-{{hemi}}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"
    output:
        func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_midthickness.func.gii"
    params:
        wb = lambda w: config["workbench_path"]
    shell:
        """
        {params.wb}/wb_command -volume-to-surface-mapping {input.volumemap} {input.surf} {output.func} -trilinear
        """

# --- Volume to Surface Mapping: white ---
rule v2s_white:
    input:
        volumemap = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/maps/{{subject}}_{{session}}_space-nativepro_map-{{feature}}.nii.gz",
        surf = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/surf/{{subject}}_{{session}}_hemi-{{hemi}}_space-nativepro_surf-fsnative_label-white.surf.gii"
    output:
        func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_white.func.gii"
    params:
        wb = lambda w: config["workbench_path"]
    shell:
        """
        {params.wb}/wb_command -volume-to-surface-mapping {input.volumemap} {input.surf} {output.func} -trilinear
        """

# --- Volume to Surface Mapping: sfwm ---
rule v2s_sfwm:
    input:
        volumemap = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/maps/{{subject}}_{{session}}_space-nativepro_map-{{feature}}.nii.gz",
        surf = "{output_dir}/{subject}/{session}/structural/{subject}_{session}_{hemi}_sfwm-{depth}mm.surf.gii"
    output:
        func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_sfwm-{depth}mm.func.gii"
    params:
        wb = lambda w: config["workbench_path"]
    shell:
        """
        {params.wb}/wb_command -volume-to-surface-mapping {input.volumemap} {input.surf} {output.func} -trilinear
        """

# --- Calculate Distances and Gradients ---
rule calc_dist_grad:
    input:
        midthickness_func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_midthickness.func.gii",
        white_func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_white.func.gii",
        sfwm_funcs = expand("{tmp_dir}/{{subject}}_{{session}}_{{hemi}}_{{feature}}_sfwm-{{depth}}mm.func.gii", depth=config['sfwm_depths']),
        midthickness_surf = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/surf/{{subject}}_{{session}}_hemi-{{hemi}}_space-nativepro_surf-fsnative_label-midthickness.surf.gii",
        white_surf = lambda w: f"{dmicapipe}/{{subject}}/{{session}}/surf/{{subject}}_{{session}}_hemi-{{hemi}}_space-nativepro_surf-fsnative_label-white.surf.gii",
        sfwm_surfs = expand("{output_dir}/{{subject}}/{{session}}/structural/{{subject}}_{{session}}_{{hemi}}_sfwm-{{depth}}mm.surf.gii", depth=config['sfwm_depths'])
    output:
        raw = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-raw.func.gii",
        dist = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-dist.func.gii",
        grad = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-grad.func.gii"
    script:
        "scripts/calc_dist_grad.py"

# --- Metric Smoothing Rule ---
rule metric_smoothing:
    input:
        surf = "{surf_file}",
        func = "{tmp_dir}/{subject}_{session}_{hemi}_{feature}_{surf_label}.func.gii"
    output:
        smooth = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_smooth-{smoothing}mm.func.gii"
    params:
        wb = lambda w: config["workbench_path"],
        smoothing = lambda w: config.get("cortical_smoothing", 5)
    shell:
        """
        {params.wb}/wb_command -metric-smoothing {input.surf} {input.func} {params.smoothing} {output.smooth}
        """

# --- Set Structure Rule ---
rule set_structure:
    input:
        smooth = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_smooth-{smoothing}mm.func.gii"
    output:
        final = "{output_dir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_smooth-{smoothing}mm_final.func.gii"
    params:
        wb = lambda w: config["workbench_path"]
    shell:
        """
        {params.wb}/wb_command -set-structure {input.smooth} CORTEX_LEFT if {{wildcards.hemi}} == "L" else CORTEX_RIGHT
        cp {input.smooth} {output.final}
        """

# Add more rules as needed for distance/gradient calculation, etc. 