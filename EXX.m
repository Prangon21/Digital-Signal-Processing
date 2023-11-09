% info = audiodevinfo;
recorder1 = audiorecorder(44100,16,1,3); 
recorder2 = audiorecorder(48000,16,1,4);

disp('Start speaking..')
record(recorder1);
record(recorder2); 
pause(5);


disp('End of Recording.');
stop(recorder1);
stop(recorder2); 
