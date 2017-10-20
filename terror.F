      subroutine terror(text)
      character * (*) text
      write(6,*)text
      stop 'error'
      end
