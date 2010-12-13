pro get_t_string

  common GLOBAL

  t_string = ' '+$
             strtrim(string(t_min),2)+$
             ' < t < '+$
             strtrim(string(t_max),2)

end
