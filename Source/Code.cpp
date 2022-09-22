#include "Code.H"
#include "Utils/WarpXUtil.H"
#include "Input/MaterialProperties.H"
#include "Input/GeometryProperties.H"



c_Code* c_Code::m_instance = nullptr;


c_Code& c_Code::GetInstance() 
{

    if (!m_instance) {
        m_instance = new c_Code();
    }
    return *m_instance;

}


void
c_Code::ResetInstance ()
{
    delete m_instance;
    m_instance = nullptr;
}


c_Code::c_Code ()
{
    m_instance = this;

    ReadData();

}


c_Code::~c_Code ()
{

 //

}


void 
c_Code::ReadData ()
{

     m_pGeometryProperties = std::make_unique<c_GeometryProperties>();
     m_pMaterialProperties = std::make_unique<c_MaterialProperties>();

}


void 
c_Code::InitData ()
{
 
    m_pGeometryProperties->InitData();
    m_pMaterialProperties->InitData();

}

