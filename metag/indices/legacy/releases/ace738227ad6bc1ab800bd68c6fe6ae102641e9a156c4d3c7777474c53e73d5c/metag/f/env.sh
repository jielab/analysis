#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./env.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# Standalone on either cluster. Select modules in the login shell or an explicit
# METAG_ENV_FILE; no TaiYi-specific compiler/Python path is assumed.
_metag_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "$_metag_root/../0f/0os.sh"
if [[ -n ${METAG_ENV_FILE:-} ]]; then
    source "$METAG_ENV_FILE"
fi
if [[ -z ${METAG_PYTHON:-} && -x /work/sph-huangj/software/venvs/metag-20261115/bin/python ]]; then
    export METAG_PYTHON=/work/sph-huangj/software/venvs/metag-20261115/bin/python
fi
if [[ -z ${METAG_PYTHON:-} && -n ${CONDA_PREFIX:-} && -x $CONDA_PREFIX/bin/python ]]; then
    export METAG_PYTHON=$CONDA_PREFIX/bin/python
fi
if [[ -z ${METAG_PYTHON:-} && -x "$_metag_root/.venv/bin/python" ]]; then
    export METAG_PYTHON="$_metag_root/.venv/bin/python"
fi
if [[ $submit == hpc ]]; then
    if [[ -z ${METAG_PYTHON:-} && -x $OS_GU_CONDA/bin/python ]]; then
        export METAG_PYTHON=$OS_GU_CONDA/bin/python
    fi
    # QiMing's system GCC 8 needs a newer C++17 filesystem runtime.
    if [[ $(g++ -dumpversion 2>/dev/null) == 8* ]]; then
        if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
            source /etc/profile.d/modules.sh
        fi
        if type module >/dev/null 2>&1; then module load gcc/12.2.0; fi
    fi
fi
if [[ -z ${METAG_PYTHON:-} && -x "$_metag_root/.venv/bin/python" ]]; then
    METAG_PYTHON="$_metag_root/.venv/bin/python"
fi
export METAG_PYTHON=${METAG_PYTHON:-python3}
export CC=${CC:-gcc} CXX=${CXX:-g++}
export OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OMP_NUM_THREADS
export PYTHONDONTWRITEBYTECODE=1

