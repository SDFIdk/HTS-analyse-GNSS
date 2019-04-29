function status = csv2txt(filename_input,filename_output)
  
  
  % Formatstring (change this if the input files get different fields.
  T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
  [REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
  textread(filename_input, T, "delimiter", ",", "headerlines", 1);
  
  % Initialize output string
  output_string = "";
  
  % 
  formatSpec = '# crt_igs08 %s \n';
  formatSpec2 = '%i %0.5f m %0.5f m %0.5f m \n';
  
  [Y, I, J] = unique (REFNR, "first");
  
    for i = 1:length(Y);
    refnr = Y(i);
    gpsnr = GPSNR(I(i));
    
    if i < length(Y)
      xkoords = XKOOR(I(i):(I(i+1)-1));
      ykoords = YKOOR(I(i):(I(i+1)-1));
      zkoords = ZKOOR(I(i):(I(i+1)-1));
      epochs = EPOCH(I(i):(I(i+1)-1));
    else
      xkoords = XKOOR(I(i):end);
      ykoords = YKOOR(I(i):end);
      zkoords = ZKOOR(I(i):end);
      epochs = EPOCH(I(i):end);
    end
  
    for j = 1:length(xkoords)
      output_string = [output_string sprintf(formatSpec,epochs{j})];
      output_string = [output_string sprintf(formatSpec2, refnr ,xkoords(j),...
                       ykoords(j),zkoords(j))];
    end    
  end
  output_string = [output_string '-1z'];
  
  fileID = fopen(filename_output,'w');
  fprintf(fileID,output_string);
  fclose(fileID);
  
  status = "done!";
end