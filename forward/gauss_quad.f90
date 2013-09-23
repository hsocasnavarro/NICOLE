Module Gauss_quad
!
  Real, Dimension(5,5) :: XMU_Gauss_Quad, WMU_Gauss_Quad
!
  Data XMU_Gauss_Quad(1,:)/0.5,0.,0.,0.,0./
  Data WMU_Gauss_Quad(1,:)/1. ,0.,0.,0.,0./
!
  Data XMU_Gauss_Quad(2,:)/0.1127017,0.7886751,0.,0.,0./
  Data WMU_Gauss_Quad(2,:)/0.5,0.5,0.,0.,0./
!
  Data XMU_Gauss_Quad(3,:)/0.1127017,0.5000000,0.8872983,0.,0./
  Data WMU_Gauss_Quad(3,:)/0.2777778,0.4444444,0.2777778,0.,0./
!
  Data XMU_Gauss_Quad(4,:)/6.9431841E-02,0.3300095,0.6699905,0.9305682,0./
  Data WMU_Gauss_Quad(4,:)/0.1739274,0.3260726,0.3260726,0.1739274,0./
!
  Data XMU_Gauss_Quad(5,:)/4.6910077E-02,0.2307653,0.5000000,.7692347,0.9530900/
  Data WMU_Gauss_Quad(5,:)/0.1184634,0.2393143,0.2844445,0.2393143,0.1184634/
!
End Module Gauss_quad
