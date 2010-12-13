PRO gyro_RZ_mpeg, data, ITMIN = itmin, ITMAX = itmax, ITSTEP = itstep, $
                  FILENAME = filename

  winsize = 256
  WINDOW, XS=winsize+16, YS=winsize

  DEFAULT, itmin, 0
  DEFAULT, itmax, 100
  DEFAULT, itstep, 2

  mpegFrame = 0
  mpegID = MPEG_OPEN([winsize,winsize], QUALITY = 75)

  FOR it = itmin, itmax, itstep DO BEGIN
      PLOT_GYRO_RZ, data, it=it, SF=128, /PLOT_N, CHARSIZE=0.9
      image = TVRD(/TRUE)
      MPEG_PUT, mpegID, IMAGE = image[*,0:winsize-1,*], FRAME = mpegFrame, $
                /ORDER
      mpegFrame += 1
  ENDFOR
  DEFAULT, filename, 'gyroRZ.mpg'
  MPEG_SAVE, mpegID, FILENAME = '/Users/cholland/Desktop/' + filename
  MPEG_CLOSE, mpegID

END ;gyro_RZ_mpeg
