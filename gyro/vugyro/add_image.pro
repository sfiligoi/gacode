pro add_image,pngfile,tag

  a=tvrd(/true)
  write_png,pngfile,a
  printf,1,'<img src='+pngfile+' alt="'+tag+'">'
  print,'Wrote '+pngfile

end
