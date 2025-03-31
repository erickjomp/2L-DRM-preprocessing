module slabvars
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, slabvars!"
  end subroutine say_hello
end module slabvars
