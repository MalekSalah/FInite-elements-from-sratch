function [P T] = lecture_solution(file_name)
  file = fopen( file_name,'r');
  if file <=0
    msg=['Le fichier de maillage : ' nomfile ' n a pas ete trouve'];
    error(msg);
  endif

  while ~strcmp(fgetl(file),'$Time'),  end
  T = str2num(fgetl(file));

  while ~strcmp(fgetl(file),'$Solution'), end
  for i=1:length(T)
    P(i,:) = str2num(fgetl(file));
  endfor
  fclose(file);
end
