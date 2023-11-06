clear all
[P T] = lecture_solution("Solution.txt");
[Nbpt,Nbtri,Coorneu,Numtri]=lecture_maillage("maillage.txt");

figure;
axis([0,2,0,2,-1,1])
trisurf(Numtri, Coorneu(:,1),Coorneu(:,2),P(1,:)); % plot initial surface
axis([0,2,0,2,0,4])
 drawnow;
shading flat
##shading faceted
##shading interp
colorbar;



for i = 1:length(T)
  trisurf(Numtri, Coorneu(:,1),Coorneu(:,2),P(i,:));
  axis([0,2,0,2,0,4])
  drawnow;                % refresh plot
  pause(0.05);             % pause to control animation speed
end


% ajouter eventuellement un titre
title(titre);
