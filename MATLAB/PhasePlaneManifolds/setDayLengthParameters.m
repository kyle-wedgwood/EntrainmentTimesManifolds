function p = setDayLengthParameters(p, day_length)
  p.day_length = day_length;
  p.length_scaling = sqrt(sin(pi/12.0*p.day_length)^2/(2.0*(1.0-cos(pi/12.0*p.day_length))));
  if (((p.day_length < 12) && (p.length_scaling > 0)) || ((p.day_length > 12) && (p.length_scaling < 0)))
       p.length_scaling = -p.length_scaling;
  end
  p.length_shift = 12.0/pi*asin(p.length_scaling);
end