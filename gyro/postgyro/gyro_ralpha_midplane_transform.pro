FUNCTION GYRO_RALPHA_MIDPLANE_TRANSFORM, field, ALPHA = alpha, $
  SCALE_FACTOR = sf
;
; C. Holland, UCSD
; v1.0: 2/21/2007- Implements existing VuGYRO algorithm
;
; takes GYRO midplane data and translates from field = complexarr(n_r,n_n) to
; field_ralpha = flatrr(n_r,n_y=sf*n_n)
;
; ALPHA is a grid vector, SF is the interpolation along alpha to use
;

  fieldsize = SIZE(field)
  n_r = fieldsize[1]
  n_n = fieldsize[2]
  IF N_ELEMENTS(sf) EQ 0 THEN sf = 1
  n_y = sf*n_n

  ;generate n->alpha phase factor for Fourier transform
  alpha = FINDGEN(n_y)/(n_y-1)
  evec = COMPLEXARR(n_n,n_y)
  cn = FLTARR(n_n) + 2.0
  cn[0] = 1.0
  FOR i_n=0,n_n-1 DO $
     evec[i_n,*] = exp(-2*!PI*COMPLEX(0,1)*i_n*alpha[*])*cn[i_n]
  
  ;do transform
  ;f(r_i,alpha_j) = Re{sum on n field(r_i,n)exp(-2*pi*n*alpha_j)}
  ;-->f(r_i,alpha_j) = Re{field(r_i,0)} + 
  ; 2*Re{field(r_i,1:n_n-1)exp(-2*pi*[1:n_n-1]*alpha_j)}
  field_ralpha = FLTARR(n_r, n_y)
  FOR ii = 0, n_r-1 DO FOR jj=0, n_y - 1 DO $
    field_ralpha[ii,jj] += TOTAL(FLOAT(field[ii,*]*evec[*,jj]))

  RETURN, field_ralpha


END ;ralpha_midplane_transform
