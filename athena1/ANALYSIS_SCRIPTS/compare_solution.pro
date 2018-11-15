pro compare_solution, start=start, stop=stop, prob=prob,dtout=dtout, xc=xc, yc=yc, zc=zc, time, time_a, density_center, density_average, density_a, radius, radius_a

if (not keyword_set(dtout)) then begin
  dtout = 0.01
endif

if (not keyword_set(nproc)) then begin
  nproc = 8
endif

if (not keyword_set(xc)) then begin
  xc = 0.0
endif

if (not keyword_set(yc)) then begin
  yc = 0.0
endif

if (not keyword_set(zc)) then begin
  zc = 0.0
endif

if (not keyword_set(prob)) then begin
  prob = "Psphere"
endif

nt = stop-start+1
time = fltarr(nt) 
density_center = fltarr(nt)
density_average = fltarr(nt)
radius = fltarr(nt)

FOR frame=start,stop DO BEGIN

  if (frame LT 10) then begin
    filename = [prob,'.000',string(frame),'.dpar.vtk']
  endif else if (frame LT 100) then begin
    filename = [prob,'.00',string(frame),'.dpar.vtk']
  endif else if (frame LT 1000) then begin
    filename = [prob,'.0',string(frame),'.dpar.vtk']
  endif else begin
    filename = [prob,'.',string(frame),'.dpar.vtk']
  endelse

  filename = strcompress(strjoin(filename),/REMOVE_ALL)

  readvtkdpar,dpar,grid,filename
  x = grid.x
  y = grid.y
  z = grid.z
  nx = n_elements(x)
  ny = n_elements(y)
  nz = n_elements(z)

  time(frame) = float(frame)*dtout
  density_center(frame) = dpar(nx/2,ny/2,nz/2)

;
; Calculate outer radius of sphere and average density
; 


; Move from left to right until there is a reasonably strong jump in the density (probably not the best way to do this...)
  for i=0,nx-1 do begin 
    r = abs(x(i)-xc)
    if (dpar(i,ny/2,nz/2) gt 0.1*density_center(frame)) then begin
       radius(frame) = r
       break
    endif
  endfor

; Now average density within this sphere
  if (frame eq start) then rp = fltarr(nx,ny,nz)
  for k=0,nz-1 do begin
  for j=0,ny-1 do begin
  for i=0,nx-1 do begin
    rp(i,j,k) = sqrt((x(i)-xc)^2.+(y(j)-yc)^2.+(z(k)-zc)^2.)
  endfor
  endfor
  endfor

  within_sphere = where(rp le radius(frame))
  density_average(frame) = AVG(dpar(within_sphere))

ENDFOR

;
; Here is the analytic solution
;

; Set 4piG to whatever it is in the code

fpg = 0.1
G = fpg/(4.*!PI)

; Set radius to whatever it is in the code
r0 = 0.75

; Set density to whatever it is in the code
d0 = 2.

; Calculate mass
M = 4./3.*!PI*r0^3.*d0

; Generate array of alpha values
alpha = findgen(1000)/1000*!PI/2.
time_a = (alpha+0.5*sin(2.*alpha))/sqrt(2.*G*M/r0^3.)
density_a = d0/cos(alpha)^6.
radius_a = r0*cos(alpha)^2.

; Calculate free fall time
tff = sqrt(3.*!PI/(32.*G*d0))

; Normalize everything
density_center /= d0
density_average /= d0
density_a /= d0
time /= tff
time_a /= tff
radius /= r0
radius_a /= r0

end
