#!/bin/bash

targSuffix=${1}

python3 IM_general_make_df.py --analysis_params IM_dediff_config.json --mode dediff --set ${targSuffix}


