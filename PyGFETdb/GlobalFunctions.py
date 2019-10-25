#
#  Author: @dragc 25/10/19 22:48.
#
#  Global Functions that do not fit in the previous files.
#


def Divide(Dividend=0, Divisor=1):
    __doc__ = "Divide (Dividend=0, Divisor=1)" \
              "" \
              "Params" \
              "------" \
              "" \
              "Dividend :" \
              "Divisor  : " \
              "" \
              "Divides Dividend and Divisor avoiding division by zero." \
              "" \
              ""
    if (Divisor == 0):
        return Dividend
    else:
        return Dividend / Divisor
