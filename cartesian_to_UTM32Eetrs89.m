function [X, Y, h] = cartesian_to_UTM32Eetrs89(X,Y,Z)
  % Transform from cartesian to UTM using proj (cct.exe)
  
  out = [X';Y';Z'];
  fileID = fopen('temp\\xyz.txt','w');
  formatSpec = '%7.5f %7.5f %7.5f\n';
  fprintf(fileID,formatSpec,out);
  fclose(fileID);

  % Using PROJ.4 / CCT.exe
  [stat,output] = system("type temp\\xyz.txt | bin\\cct.exe +proj=pipeline +ellps=GRS80 +step +proj=cart +inv +step +proj=utm +zone=32");

  C = strsplit(output);
  
  X = [];
  Y = [];
  h = [];
  
  for i = 1:4:length(C)-2
    X = [X ; str2num(C{i})];
    Y = [Y ; str2num(C{i+1})];
    h = [h ; str2num(C{i+2})];
  
end