#ifndef TEST_MATH_LIBRARY_H_
#define TEST_MATH_LIBRARY_H_

#include "Tests/All_Tests.H"

template <typename Test=Test_MM_Mul>
class Test_MathLibrary
{
   private:
    Test m_test;

   public:
    Test_MathLibrary(Test&& test_instance) : m_test(std::move(test_instance))
    {
        m_test.Define();
        m_test.Initialize();
        m_test.Print_Input();
    }

    void Perform_Test()
    {
        m_test.Perform_Test();
    }

    void Print_Output() {
        m_test.Print_Output();       
    }

};

#endif
