module patmo_sparsity
contains

  !compute IA and JA for sparsity
  subroutine computeSparsity()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::Ms(neqAll,neqAll),i,j,k
    integer::nnz,nnzi,jatmp(neqAll**2)

    Ms(:,:) = 0

    !chemistry sparsity
    do i=1,cellsNumber
      do j=1,reactionsNumber
        Ms((indexReactants1(j)-1)*cellsNumber+i, &
            (indexReactants1(j)-1)*cellsNumber+i) = 1
        Ms((indexReactants2(j)-1)*cellsNumber+i, &
            (indexReactants1(j)-1)*cellsNumber+i) = 1

        Ms((indexProducts1(j)-1)*cellsNumber+1, &
            (indexReactants1(j)-1)*cellsNumber+i) = 1
        Ms((indexProducts2(j)-1)*cellsNumber+1, &
            (indexReactants1(j)-1)*cellsNumber+i) = 1
        Ms((indexProducts3(j)-1)*cellsNumber+1, &
            (indexReactants1(j)-1)*cellsNumber+i) = 1

        Ms((indexReactants1(j)-1)*cellsNumber+i, &
            (indexReactants2(j)-1)*cellsNumber+i) = 1
        Ms((indexReactants2(j)-1)*cellsNumber+i, &
            (indexReactants2(j)-1)*cellsNumber+i) = 1

        Ms((indexProducts1(j)-1)*cellsNumber+1, &
            (indexReactants2(j)-1)*cellsNumber+i) = 1
        Ms((indexProducts2(j)-1)*cellsNumber+1, &
            (indexReactants2(j)-1)*cellsNumber+i) = 1
        Ms((indexProducts3(j)-1)*cellsNumber+1, &
            (indexReactants2(j)-1)*cellsNumber+i) = 1

      end do
    end do

    !diffusion sparsity
    !dni(j)/dni(j+1)/dt
    do i=1,speciesNumber
      do j=1,cellsNumber-1
        Ms((i-1)*cellsNumber+j,(i-1)*cellsNumber+j+1) = 1
      end do
    end do

    !diffusion sparsity
    !dni(j)/dni(j-1)/dt
    do i=1,speciesNumber
      do j=2,cellsNumber
        Ms((i-1)*cellsNumber+j,(i-1)*cellsNumber+j-1) = 1
      end do
    end do

    !Tgas does not change
    Ms(:,positionTgas) = 0
    Ms(positionTgas,:) = 0

    !dummy does not change
    Ms(:,positionDummy) = 0
    Ms(positionDummy,:) = 0

    !diagonal is dense by construction
    do i=1,neqAll
      Ms(i,i) = 1
    end do

    !convert the sparsity map to Yale sparisty format
    iaSparsity(1) = 1 !first row position index is 1 by definition
    nnz = 0 !count non-zero elements
    !loop on ALL the columns
    do j=1,neqAll
      nnzi = 0 !count the non zero-elements in the jth column
      !loop on the rows to find non-zero elements
      do k=1,neqAll
        !check if the element is non-zero
        if(Ms(k,j)==1) then
          nnz = nnz + 1 !increase the total non-zero elements
          nnzi = nnzi + 1 !increase the column non-zero elements
          !position of the non-zero element in the column
          jatmp(nnz) = k
        end if
      end do
      !set the integrated number of non-zero elements
      ! for the jth column (note that ia(1)=1, so j+1)
      iaSparsity(j+1) = iaSparsity(j) + nnzi
    end do

    !allocate jaa depending on the non-zero elements found
    allocate(jaSparsity(nnz))

    !copy the non-zero positions from the temporary array
    ! into the allocated array
    jaSparsity(:) = jatmp(1:nnz)

    !set non-zero elements
    nonZeroElements = nnz

  end subroutine computeSparsity

end module patmo_sparsity
