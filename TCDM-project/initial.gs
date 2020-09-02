'reinit'
'set display color white'
'c'
'open initial.ctl'
'set lon 130 140'
'set lat 20 33'

n = 1
while(n<14)
'set cmin 3600'
'set cint 25'
'set t 'n''
*'define h = hH-hM'
n = n+1
'd hhM'
pull pause
'c'
endwhile

