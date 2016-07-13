subroutine interpgene(gene,mask,nx,ny,nz,l,genei)

! This subroutine is a nearest-neighbour type interpolator. It only interpolate when
! l% (l fraction) of the nearest-neighbours express the valid value.
! used in the context of gene expression


   implicit none
   integer, intent(in) :: nx,ny,nz                   !dimensions of gene expression
   real*8, intent(in), dimension(nx,ny,nz) :: gene   !gene expression
   logical, intent(in), dimension(nx,ny,nz) :: mask  !mask for gene expression
   real*8,intent(in) :: l                            !interpolate if l% of neighbours express gene
   real*8, dimension(nx,ny,nz) :: genei              !interpolated gene expression
   integer :: i,j,k,i2,j2,k2                         !variables for loops
   real*8 :: meangene
   integer :: goodneigh,maskneigh


   genei=gene

   !iterate through elements
   do i=1,nx
   do j=1,ny
   do k=1,nz
     !skip voxels outside mask and non-negative gene expression
     if (.not.mask(i,j,k)) then
       cycle
     end if
     if (.not.(gene(i,j,k)<0)) then
       cycle
     end if
     
     goodneigh=0 ; maskneigh=0 ; meangene=0.0
     !Iterate through nearest neighbours
     do i2=(i-2),(i+2)
        !skip elements outside array index
        if (i2<1) then
          cycle
        end if
        if (i2>nx) then
          cycle
        end if
     do j2=(j-2),(j+2)
        !skip elements outside array index
        if (j2<1) then
          cycle
        end if
        if (j2>ny) then
          cycle
        end if
     do k2=(k-2),(k+2)
        !skip elements outside array index
        if (k2<1) then
          cycle
        end if
        if (k2>nz) then
          cycle
        end if
        
        !If neighbour index is outside mask, skip
        if (.not.mask(i2,j2,k2)) then
          cycle
        end if

        !count number of neighbours in mask
        maskneigh=maskneigh+1
        
        !count number of neighbours with gene expression
        !take these neighbours for mean
        if (.not.(gene(i2,j2,k2)<0)) then
          goodneigh=goodneigh+1
          meangene=meangene+gene(i2,j2,k2)
        end if
         
     end do
     end do
     end do

     !If l% of the neighbours in mask express gene
     !Fill voxel with average of its expressing neighbours
     if (REAL(goodneigh)>(REAL(maskneigh)*l)) then
         genei(i,j,k)=meangene/REAL(goodneigh)
     end if
   end do
   end do
   end do


end subroutine interpgene
