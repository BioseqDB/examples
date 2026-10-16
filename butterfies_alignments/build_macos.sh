#!/usr/bin/env bash
# Tested toolchain: Apple Silicon, macOS, LLVM 22, PostgreSQL 17.
# Prerequisites: brew install cmake ninja llvm@22 postgresql@17 pkg-config gettext
# Also bootstrap vcpkg and export VCPKG_ROOT to its checkout.
# This builds locally; it does not replace a running server or install files.
set -euo pipefail

if [[ "$(uname -s)" != Darwin || "$(uname -m)" != arm64 ]]; then
  echo "This helper describes the tested Apple Silicon toolchain. Use your platform's BioseqDB source-build instructions with BIOSEQ_BUILD_EXTERNAL_API_CLIENTS=ON." >&2
  exit 1
fi
: "${BIOSEQ_SOURCE:?Set BIOSEQ_SOURCE to your patched BioseqDB source checkout}"
: "${VCPKG_ROOT:?Set VCPKG_ROOT to a bootstrapped vcpkg checkout}"
BIOSEQ_TUTORIAL_BUILD="${BIOSEQ_TUTORIAL_BUILD:-$BIOSEQ_SOURCE/build/butterflies-tutorial}"
BIOSEQ_TUTORIAL_LLVM="$(brew --prefix llvm@22)"
BIOSEQ_TUTORIAL_PG="$(brew --prefix postgresql@17)"

cmake -S "$BIOSEQ_SOURCE" -B "$BIOSEQ_TUTORIAL_BUILD" -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  -DCMAKE_C_COMPILER="$BIOSEQ_TUTORIAL_LLVM/bin/clang" \
  -DCMAKE_CXX_COMPILER="$BIOSEQ_TUTORIAL_LLVM/bin/clang++" \
  -DCMAKE_BUILD_TYPE=Debug \
  '-DCMAKE_C_FLAGS=-O3 -DNDEBUG -fno-omit-frame-pointer' \
  '-DCMAKE_CXX_FLAGS=-O3 -DNDEBUG -fno-omit-frame-pointer' \
  -DCMAKE_C_FLAGS_DEBUG= -DCMAKE_CXX_FLAGS_DEBUG= \
  -DCMAKE_OSX_ARCHITECTURES=arm64 -DVCPKG_TARGET_TRIPLET=arm64-osx \
  -DPostgreSQL_CONFIG="$BIOSEQ_TUTORIAL_PG/bin/pg_config" \
  -DBIOSEQ_BUILD_EXTERNAL_API_CLIENTS=ON \
  -DBIOSEQ_BUILD_WEBSITE=OFF -DBIOSEQ_BUILD_PAPERS=OFF \
  -DBIOSEQ_BUILD_CACHE=ON -DBIOSEQ_TEST_INDEXES=OFF \
  -DBIOSEQ_TEMPLATE_FLATTEN_PROFILE=all_alphabets_uint64 \
  -DBIOSEQ_HEAVY_JOBS=2 \
  -DBIOSEQ_CODEGEN_CXX_COMPILER="$BIOSEQ_TUTORIAL_LLVM/bin/clang++" \
  -DBIOSEQ_CODEGEN_C_COMPILER="$BIOSEQ_TUTORIAL_LLVM/bin/clang" \
  -DBIOSEQ_CODEGEN_CLANG_BINARY="$BIOSEQ_TUTORIAL_LLVM/bin/clang++" \
  -DBIOSEQ_CODEGEN_LLVM_CONFIG="$BIOSEQ_TUTORIAL_LLVM/bin/llvm-config" \
  -DBIOSEQ_CODEGEN_CLANG_SCAN_DEPS="$BIOSEQ_TUTORIAL_LLVM/bin/clang-scan-deps" \
  -DBIOSEQ_CODEGEN_STDLIB_MODULES_JSON="$BIOSEQ_TUTORIAL_LLVM/lib/c++/libc++.modules.json" \
  -DBIOSEQ_CODEGEN_CLANG_STDLIB_MODULES_JSON="$BIOSEQ_TUTORIAL_LLVM/lib/c++/libc++.modules.json"
cmake --build "$BIOSEQ_TUTORIAL_BUILD" \
  --target postgresql_bioseq postgresql_bioseq_cache --parallel 3
echo "Built in $BIOSEQ_TUTORIAL_BUILD. Install into a development PostgreSQL instance with:"
echo "cmake --install \"$BIOSEQ_TUTORIAL_BUILD\" --component bioseq"
