class zbenv:
    def __init__(self, connectome_workbench_path, num_threads=None):
        """
        Initialize the z-brains environment configuration.
        
        Parameters:
        -----------
        connectome_workbench_path : str
            Path to the Connectome Workbench executables
        num_threads : int or None, default=None
            Number of threads to use for Connectome Workbench operations.
            If None, uses system default (typically all available cores).
        """
        self.connectome_workbench_path = connectome_workbench_path
        self.num_threads = num_threads
        
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
            env['OMP_NUM_THREADS'] = str(self.num_threads)
            
        return env