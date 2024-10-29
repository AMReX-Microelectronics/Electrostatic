# For clang-format
find ../../Source -type f \( -name "*.cpp" -o -name "*.H" \) -exec clang-format -i {} \;

# For sed
find ../../Source -type f \( -name "*.cpp" -o -name "*.H" \) -exec sed -i 's/[[:blank:]]\+$//' {} \;
