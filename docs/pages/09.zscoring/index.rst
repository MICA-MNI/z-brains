.. _normz:

.. title:: Z-scoring

Z-scoring
============================================================

This section describes how to run the z-scoring module to normalize single patient data against a cohort of controls.


Setting up your control dataset
--------------------------------------------------------
You will need to set up a csv file with two mandatory columns: 
One column containing the subject IDs for all controls to be included in the analysis (column name: ``ID``), 
as well as one column with corresponding session name to be analyzed (column name: ``SES``). 
Note: if your dataset does not contain session information, then leave the ``SES`` column blank. 

An example control table can be found in our GitHub repository.

Regional analysis
--------------------------------------------------------

After defining your directories as in the previous section (``-proc``), you can run the regional z-scoring analysis as follows:

    .. code-block:: bash
        :caption: Basic z-brains run: hippocampal processing
        :linenos:
        
        # Define path to control dataset csv
        hc_info=/PATH/TO/CSV/control.csv
        
        # Information on patient to be analyzed
        px_id="PX001"
        px_ses="01"

        z-brains -sub "$id" -ses "$ses" \
            -rawdir "${rawdir}" \
            -micapipedir "${micapipedir}" \
            -hippdir "${hippdir}" \
            -outdir "${outdir}" \
            -run regional \
            -approach "zscore" \
            -demo_cn "${hc_info}" \
            -verbose 2

Asymmetry analysis
--------------------------------------------------------

Similar to the ``regional`` analysis, you can generated z-scored asymmetry maps with the following approach:

    .. code-block:: bash
        :caption: Basic z-brains run: hippocampal processing
        :linenos:
        
        # Define path to control dataset csv
        hc_info=/PATH/TO/CSV/control.csv
        
        # Information on patient to be analyzed
        px_id="PX001"
        px_ses="01"

        z-brains -sub "$id" -ses "$ses" \
            -rawdir "${rawdir}" \
            -micapipedir "${micapipedir}" \
            -hippdir "${hippdir}" \
            -outdir "${outdir}" \
            -run asymmetry \
            -approach "zscore" \
            -demo_cn "${hc_info}" \
            -verbose 2
