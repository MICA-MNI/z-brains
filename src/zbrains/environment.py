class zbenv:
    def __init__(
        self,
        connectome_workbench_path,
        num_threads=1,
        num_threads_wb=1,
        normalization_jobs=None,
        normalization_threads_per_job=None,
        ravel_tmp_dir=None,
    ):
        """
        Initialize the z-brains environment configuration.
        
        Parameters:
        -----------
        connectome_workbench_path : str
            Path to the Connectome Workbench executables
        num_threads : int or None, default=None
            Total CPU budget used by the processing pipeline.
        normalization_jobs : int or None, default=None
            Concurrent subject-level normalization jobs. When omitted, z-brains
            selects a bounded value from ``num_threads``.
        normalization_threads_per_job : int or None, default=None
            ANTs/SynthSeg threads assigned to each normalization job. When
            omitted, the total CPU budget is divided between concurrent jobs.
        ravel_tmp_dir : path-like or None, default=None
            Local scratch directory for transient RAVEL memmaps. The
            ``ZBRAINS_RAVEL_TMPDIR`` or ``SLURM_TMPDIR`` environment variable is
            used when this is omitted.
        """
        self.connectome_workbench_path = connectome_workbench_path
        self.num_threads = num_threads
        self.num_threads_wb = num_threads_wb
        self.normalization_jobs = normalization_jobs
        self.normalization_threads_per_job = normalization_threads_per_job
        self.ravel_tmp_dir = ravel_tmp_dir
        
    def get_workbench_env(self):
        """
        Get environment variables for Connectome Workbench operations.
        
        Returns:
        --------
        dict
            Environment variables to use with subprocess calls
        """
        import os
        env = os.environ.copy()
        
        if self.num_threads is not None:
            env['OMP_NUM_THREADS'] = str(self.num_threads_wb)
            
        return env
