function [Nbpt,Nbtri,Coorneu,Numtri]=lecture_maillage(file_name)

  file = fopen( file_name,'r');
  if file <=0
    msg=['Le fichier de maillage : ' nomfile ' n a pas ete trouve'];
    error(msg);
  endif

  while ~strcmp(fgetl(file),'$Nodes'),  end
  Nbpt = str2num(fgetl(file));
  Coorneu = zeros(Nbpt,2);

  for i=1:Nbpt
    tmp= str2num(fgetl(file));
    Coorneu(i,:) = tmp;
  endfor

  while ~strcmp(fgetl(file),'$Elements'), end
  Nbtri = str2num(fgetl(file));
  Numtri = zeros(Nbtri,3);

  for i=1:Nbtri
    tmp= str2num(fgetl(file));
    Numtri(i,:) = tmp;
  endfor
  fclose(file);
end
