#ifndef TEST_MATH_LIBRARY_H_
#define TEST_MATH_LIBRARY_H_

#include <AMReX_Print.H>

#include "Tests/All_Tests.H"

template <typename Test = Test_MM_Mul>
class Test_MathLibrary
{
   private:
    Test m_test;
    int repetitions = 1;
    std::vector<int> print_flags;

   public:
    Test_MathLibrary(Test&& test_instance, const std::vector<int>& flags)
        : m_test(std::move(test_instance)), print_flags(flags)
    {
        m_test.Define();
        m_test.Initialize();

        if (print_flags[0]) m_test.Print_Input();
    }

    void Perform_Test(int repetitions = 1)
    {
        for (int i = 0; i < repetitions; ++i)
        {
            m_test.Perform_Test();
        }

        if (print_flags[1]) m_test.Print_Output();

        m_test.Verify(print_flags[2]);

        m_test.Clear();
    }
};

#endif
