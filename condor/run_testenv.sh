#!/bin/bash

# ----------------------------
source /storage/fdunphome/wangchunzheng/miniconda3/etc/profile.d/conda.sh
conda activate cpp_dev

export CC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc"
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++"
unset ROOTSYS

export PATH=$(echo ":$PATH:" | sed -E 's#:/opt/root61404/bin:#:#g' | sed 's#^:##;s#:$##')
export LD_LIBRARY_PATH=$(echo ":$LD_LIBRARY_PATH:" | sed -E 's#:/opt/root61404/lib:#:#g' | sed 's#^:##;s#:$##')

export PATH="$CONDA_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="/storage/fdunphome/wangchunzheng/CWCoalProject/lib64:$LD_LIBRARY_PATH"


echo "=== ENV CHECK ==="
echo "which cc: $(which $CC)"
echo "which c++: $(which $CXX)"
echo "conda env  : $CONDA_DEFAULT_ENV"
echo "PATH       : $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

echo "=== RUNNING TEST PROGRAM ==="
echo '#include<iostream>' > test.cpp
echo 'int main() { std::cout << "Hello from Condor!" << std::endl; return 0; }' >> test.cpp
$CXX test.cpp -o test.out && ./test.out || echo "Compile or run failed"


conda deactivate
echo "=== JOB COMPLETED ==="
