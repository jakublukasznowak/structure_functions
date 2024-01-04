function d = mydatapath

comptype=computer;
if strcmp(comptype,'PCWIN64')
    main='C:\jnowak';
elseif strcmp(comptype,'GLNXA64')
    main='/home/pracownicy/jnowak';
else
    error('Operating system unknown')
end

d = [main,filesep,'NEXTGEMS',filesep,'data',filesep,'ATR'];

end