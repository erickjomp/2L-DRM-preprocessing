module IO
    implicit none
contains
    ! if not removing initial or last file, rem_file should equal to 0
    subroutine read_list_files(text_file, files, rem_file) 
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/4
        ! type :: vlstring_t
        ! character(len=:), allocatable :: string
        ! end type
        integer :: iun, nrecs, istat
        character(len=*), intent(in)  :: text_file
        character(len=200), intent(out), allocatable :: files(:)
        integer, intent(in)               :: rem_file
        character(len=200) :: buffer
        character(len=1) :: dummy
        integer          :: i
        character(len=200), allocatable :: files_temp(:)

        open(newunit=iun, file=text_file, form="formatted", status = "old")
        ! count number of records (note you don't have to read the line of text just the first char)
        nrecs = 0
        count_loop: do
            read(iun, '(a1)', iostat=istat) dummy
            if (istat /= 0) exit count_loop
            nrecs = nrecs+1
        end do count_loop

        allocate(files(nrecs))
        rewind(iun)
        read_loop: do i=1,nrecs
            buffer=repeat(" ", 200)
            read(iun,'(a)', iostat=istat) buffer
            ! print *,buffer
            if (istat /= 0) then
                print *,' error while reading text file'
                exit read_loop
            end if
            files(i) = trim(buffer)
        end do read_loop
        close(iun)

        ! only in case removing file
        if (rem_file == 1 .or. rem_file == 2) then
            allocate(files_temp(nrecs - 1))
            if (rem_file == 1) then
                files_temp = files(2:)
            else if (rem_file == 2) then
                files_temp = files(:size(files)-1)
            end if
            deallocate(files)
            allocate(files(nrecs-1))
            files = files_temp
            deallocate(files_temp)
        end if
    end subroutine

    ! https://fortran-lang.github.io/fpm/proc/basename.html
    function get_basename(path,suffix) result (base)

        character(*), intent(In) :: path
        logical, intent(in), optional :: suffix
        character(:), allocatable :: base
    
        character(:), allocatable :: file_parts(:)
        logical :: with_suffix
    
        if (.not.present(suffix)) then
            with_suffix = .true.
        else
            with_suffix = suffix
        end if
    
        call split(path,file_parts,delimiters='\/')
        if(size(file_parts)>0)then
           base = trim(file_parts(size(file_parts)))
        else
           base = ''
        endif
        if(.not.with_suffix)then
            call split(base,file_parts,delimiters='.')
            if(size(file_parts)>=2)then
               base = trim(file_parts(size(file_parts)-1))
            endif
        endif
    
    end function get_basename


    ! https://fortran-lang.github.io/fpm/proc/split.html#src
    subroutine split(input_line,array,delimiters,order,nulls)
        !! given a line of structure " par1 par2 par3 ... parn " store each par(n) into a separate variable in array.
        !!
        !! * by default adjacent delimiters in the input string do not create an empty string in the output array
        !! * no quoting of delimiters is supported
        character(len=*),intent(in)              :: input_line  !! input string to tokenize
        character(len=*),optional,intent(in)     :: delimiters  !! list of delimiter characters
        character(len=*),optional,intent(in)     :: order       !! order of output array sequential|[reverse|right]
        character(len=*),optional,intent(in)     :: nulls       !! return strings composed of delimiters or not ignore|return|ignoreend
        character(len=:),allocatable,intent(out) :: array(:)    !! output array of tokens
    
        integer                       :: n                      ! max number of strings INPUT_LINE could split into if all delimiter
        integer,allocatable           :: ibegin(:)              ! positions in input string where tokens start
        integer,allocatable           :: iterm(:)               ! positions in input string where tokens end
        character(len=:),allocatable  :: dlim                   ! string containing delimiter characters
        character(len=:),allocatable  :: ordr                   ! string containing order keyword
        character(len=:),allocatable  :: nlls                   ! string containing nulls keyword
        integer                       :: ii,iiii                ! loop parameters used to control print order
        integer                       :: icount                 ! number of tokens found
        integer                       :: ilen                   ! length of input string with trailing spaces trimmed
        integer                       :: i10,i20,i30            ! loop counters
        integer                       :: icol                   ! pointer into input string as it is being parsed
        integer                       :: idlim                  ! number of delimiter characters
        integer                       :: ifound                 ! where next delimiter character is found in remaining input string data
        integer                       :: inotnull               ! count strings not composed of delimiters
        integer                       :: ireturn                ! number of tokens returned
        integer                       :: imax                   ! length of longest token
    
        ! decide on value for optional DELIMITERS parameter
        if (present(delimiters)) then                                     ! optional delimiter list was present
            if(delimiters/='')then                                       ! if DELIMITERS was specified and not null use it
                dlim=delimiters
            else                                                           ! DELIMITERS was specified on call as empty string
                dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0) ! use default delimiter when not specified
            endif
        else                                                              ! no delimiter value was specified
            dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)    ! use default delimiter when not specified
        endif
        idlim=len(dlim)                                                   ! dlim a lot of blanks on some machines if dlim is a big string
    
        if(present(order))then; ordr=lower(adjustl(order)); else; ordr='sequential'; endif ! decide on value for optional ORDER parameter
        if(present(nulls))then; nlls=lower(adjustl(nulls)); else; nlls='ignore'    ; endif ! optional parameter
    
        n=len(input_line)+1                        ! max number of strings INPUT_LINE could split into if all delimiter
        allocate(ibegin(n))                        ! allocate enough space to hold starting location of tokens if string all tokens
        allocate(iterm(n))                         ! allocate enough space to hold ending location of tokens if string all tokens
        ibegin(:)=1
        iterm(:)=1
    
        ilen=len(input_line)                                           ! ILEN is the column position of the last non-blank character
        icount=0                                                       ! how many tokens found
        inotnull=0                                                     ! how many tokens found not composed of delimiters
        imax=0                                                         ! length of longest token found
    
        select case (ilen)
    
        case (0)                                                      ! command was totally blank
    
        case default                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
            icol=1                                                      ! initialize pointer into input line
            INFINITE: do i30=1,ilen,1                                   ! store into each array element
                ibegin(i30)=icol                                         ! assume start new token on the character
                if(index(dlim(1:idlim),input_line(icol:icol))==0)then  ! if current character is not a delimiter
                iterm(i30)=ilen                                       ! initially assume no more tokens
                do i10=1,idlim                                        ! search for next delimiter
                    ifound=index(input_line(ibegin(i30):ilen),dlim(i10:i10))
                    IF(ifound>0)then
                        iterm(i30)=min(iterm(i30),ifound+ibegin(i30)-2)
                    endif
                enddo
                icol=iterm(i30)+2                                     ! next place to look as found end of this token
                inotnull=inotnull+1                                   ! increment count of number of tokens not composed of delimiters
                else                                                     ! character is a delimiter for a null string
                iterm(i30)=icol-1                                     ! record assumed end of string. Will be less than beginning
                icol=icol+1                                           ! advance pointer into input string
                endif
                imax=max(imax,iterm(i30)-ibegin(i30)+1)
                icount=i30                                               ! increment count of number of tokens found
                if(icol>ilen)then                                     ! no text left
                exit INFINITE
                endif
            enddo INFINITE
    
        end select
    
        select case (trim(adjustl(nlls)))
        case ('ignore','','ignoreend')
            ireturn=inotnull
        case default
            ireturn=icount
        end select
        allocate(character(len=imax) :: array(ireturn))                ! allocate the array to return
        !allocate(array(ireturn))                                       ! allocate the array to turn
    
        select case (trim(adjustl(ordr)))                              ! decide which order to store tokens
        case ('reverse','right') ; ii=ireturn ; iiii=-1                ! last to first
        case default             ; ii=1       ; iiii=1                 ! first to last
        end select
    
        do i20=1,icount                                                ! fill the array with the tokens that were found
            if(iterm(i20)<ibegin(i20))then
                select case (trim(adjustl(nlls)))
                case ('ignore','','ignoreend')
                case default
                array(ii)=' '
                ii=ii+iiii
                end select
            else
                array(ii)=input_line(ibegin(i20):iterm(i20))
                ii=ii+iiii
            endif
        enddo
    end subroutine split


    ! https://fortran-lang.github.io/fpm/proc/lower.html
    elemental pure function lower(str,begin,end) result (string)

        character(*), intent(In)     :: str
        character(len(str))          :: string
        integer,intent(in),optional  :: begin, end
        integer                      :: i
        integer                      :: ibegin, iend
        string = str

        ibegin = 1
        if (present(begin))then
            ibegin = max(ibegin,begin)
        endif

        iend = len_trim(str)
        if (present(end))then
            iend= min(iend,end)
        endif

        do i = ibegin, iend                               ! step thru each letter in the string in specified range
            select case (str(i:i))
            case ('A':'Z')
                string(i:i) = char(iachar(str(i:i))+32)     ! change letter to miniscule
            case default
            end select
        end do

    end function lower

end module IO

