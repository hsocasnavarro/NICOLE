! This routine propagates forward the inputs through a rectangular Artificial
! Neural Network, and returns the result at the output nodes.
! W(l, i, j) represents the synaptic strength from neuron j in the
! layer l-1 to neuron i in layer l. Beta(l, i) is the bias added to the signal
! in neuron i, layer l. nonlin(l) is an integer vector whose elements are 0 
! if l is a linear layer, or 1 if l is a non-linear layer (tanh is then used
! as activation function). input and output are the inputs and outputs vectors,
! of size ninputs and noutputs, respectively. nperlayer is the number of
! neurons per layer and nlayers is the number of layers (note: the input
! is _NOT_ considered a layer). If the keyword Y is present, the
! neuron values are returned.
!
! Note that nmaxperlayer _MUST BE_ larger than both ninputs and noutputs!!
!
Subroutine ANN_forward(W, Beta, Nonlin, Input, Output, nlayers, &
     nmaxperlayer, nperlayer, ninputs, noutputs, y)
Implicit None
Real (Kind=8), Parameter :: a=1.7159, b=0.666666, bovera=0.388523, asq=2.94431
!Real (Kind=8), Parameter :: a=1., b=1., bovera=1., asq=1.
Integer :: nlayers, nmaxperlayer, ninputs, noutputs, i, j, k, l
Integer, Dimension(nlayers) :: Nonlin
Integer, Dimension(0:nlayers) :: nperlayer
Real (Kind=8), Dimension(nlayers, nmaxperlayer, nmaxperlayer) :: W
Real (Kind=8), Dimension(0:nlayers, nmaxperlayer) :: y
Real (Kind=8), Dimension(nlayers, nmaxperlayer) :: Beta
Real (Kind=8), Dimension(ninputs) :: Input
Real (Kind=8), Dimension(noutputs) :: Output
!
! Set input values
!
y(0, 1:ninputs)=Input(1:ninputs)
!
! Propagate forward
!
Do l=1, nlayers
!
! Do the matrix multiplications using the intrinsic F90 function
!
!   y(l, 1:nperlayer(l))=Matmul(W(l, 1:nperlayer(l), 1:nperlayer(l-1)), &
!        y(l-1, 1:nperlayer(l-1)))+Beta(l, 1:nperlayer(l))
!
! Spell out the matrix multiplications (some compilers have trouble
!    optimizing this function)
!
   Do j=1, nperlayer(l)
      y(l, j)=0.
      Do k=1, nperlayer(l-1)
         y(l, j)=y(l, j)+W(l, j, k)*y(l-1, k)
      End do
      y(l, j)=y(l, j)+Beta(l, j)
   End do
!
   If (Nonlin(l) .ne. 0) & ! It's a non-linear layer
        y(l, 1:nperlayer(l))=a*Tanh(b*y(l, 1:nperlayer(l)))
End do
!
! Network output
!
Output(1:noutputs)=y(nlayers, 1:noutputs)
!
! Done
!
Return
!
End Subroutine ANN_forward
