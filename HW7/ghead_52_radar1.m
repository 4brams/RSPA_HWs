
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subprogram name     : headr.m                               %
        % Date & time         : Mar. 8 1996 20:15                     %
        % Version             : 1                                     %
        % Purpose             : header read from binary raw data file %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fid=fopen(infile,'r');
%         fod=fopen(outfile,'w');
%         fod1=fopen(outfile2,'w');
        fod=fopen(infile,'r');
        fod1=fopen(infile,'r');
        s0=fgetl(fid);           %#GS UD3 Raw Data file
        s1=fgetl(fid);           %DEBUG plugin
        s2=fgetl(fid);           %ACQ_TIME
        s3=fgetl(fid);           %ANTENNA_COORDS
        s3_0=fgetl(fid);           %BEAM_DIRECTION
        s3_1=fgetl(fid);           %BEAM_INDEX
        s3_2=fgetl(fid);           %BEAM_LIST
        s4=fgetl(fid);           %CABLE_PHASES
%        s5=fgetl(fid);           %CAMPAIGN
        s6=fgetl(fid);           %CHANNELS
        s7=fgetl(fid);           %COMMENT
        s7_1=fgetl(fid);           %DIG_PEAK_VALUE
        s7_2=fgetl(fid);           %DIG_WIDTH
        s8=fgetl(fid);           %EXPTAG
        s9=fgetl(fid);           %FDI_OFFSETS
        s10=fgetl(fid);         %FREQUENCY
        s11=fgetl(fid);         %GATES
        s12=fgetl(fid);         %GPS_STATUS
        s13=fgetl(fid);         %INTEGRATIONS
        s14=fgetl(fid);         %INTERVAL
        s14_1=fgetl(fid);         %IPLO_LEADS_QPLO
        s15=fgetl(fid);         %LO_FREQUENCY
        s16=fgetl(fid);         %LOCATION
        s17=fgetl(fid);         %PHASE_OFFSETS
        s18=fgetl(fid);         %PLUGIN_PATH
        s19=fgetl(fid);         %PRF
        s20=fgetl(fid);         %PULSE_AMPLITUDE
        s21=fgetl(fid);         %PULSE_CODE
        s22=fgetl(fid);         %PULSE_LENGTH
        s23=fgetl(fid);         %PULSE_NUMBER
        s23_1=fgetl(fid);         %PULSE_PHASE
        s24=fgetl(fid);         %PULSE_RANGEZERO
        s25=fgetl(fid);         %PULSE_RXFILTER
        s26=fgetl(fid);         %PULSE_RXGAIN
        s27=fgetl(fid);         %PULSE_SHAPE
        s28=fgetl(fid);         %RANGE
        s29=fgetl(fid);         %RANGE_ZERO
        s30=fgetl(fid);         %RDAS_ID
        s31=fgetl(fid);         %RECEIVER_PHASES
        s32=fgetl(fid);         %RESOLUTION
        s33=fgetl(fid);         %RESULTS
        s34=fgetl(fid);         %RX_FILTER
        s35=fgetl(fid);         %RX_GAIN
%        s36=fgetl(fid);         %RX_MAXGAIN
%        s37=fgetl(fid);         %RX_MINGAIN
%        s38=fgetl(fid);         %SEQUENCE
        s39=fgetl(fid);         %SITENAME
        s40=fgetl(fid);         %TIME_ACCURACY
        s41=fgetl(fid);         %TX_PHASEFLIP
        s42=fgetl(fid);         %CONFIG
        s43=fgetl(fid);         %DATA
        
        disp(s1);
        disp(s2);
        disp(s3);
        disp(s3_0);
        disp(s4);
 %       disp(s5);
        disp(s6);
        disp(s7);
        disp(s8);
        disp(s9);
        disp(s10);
        disp(s11);
        disp(s12);
        disp(s13);
        disp(s14);
        disp(s15);
        disp(s16);
        disp(s17);
        disp(s18);
        disp(s19);
        disp(s20);
        disp(s21);
        disp(s22);
        disp(s23);
        disp(s24);
        disp(s25);
        disp(s26);
        disp(s27);
        disp(s28);
        disp(s29);
        disp(s30);
        disp(s31);
        disp(s32);
        disp(s33);
        disp(s34);
        disp(s35);
%        disp(s36);
%        disp(s37);
%        disp(s38);
        disp(s39);
        disp(s40);
        disp(s41);
        disp(s42);
        disp(s43);



%        date1=[s1(32),s1(33),s1(34),s1(35),s1(36),s1(37),s1(38),s1(39)];
        time1=str2num(s43(6:15));
        unix_time_start=datenum('1-Jan-1970');        %  days of 0000/00/00~1969/12/31 %
        days=time1/60/60/24+unix_time_start;        % days of 0000/00/00~data time %
        year_date_time=datestr(days,31);       %  date of data %
        d_time=str2num(year_date_time(12:13))+str2num(year_date_time(15:16))/60+str2num(year_date_time(18:19))/60/60;        % 24 hours %
        mon=[year_date_time(6:7)];
        day=[year_date_time(9:10)];
        year=[year_date_time(1:4)];
        mon0=str2num(mon);
        day0=str2num(day);
        year0=str2num(year);
%        date1=[s1(32),s1(33),s1(34),s1(35),s1(36),s1(37),num2str(year0)];
        sh0=year_date_time(12:13);
        sm0=year_date_time(15:16);
        ss0=year_date_time(18:19);
        h0=str2num(sh0);
        m0=str2num(sm0);
        s0=str2num(ss0);

        ls3_0=length(s3_0);
        ls3_1=length(s3_1);
        
        ls6=length(s6);
        nch=str2num(s6(9:ls6));
        
        ls19=length(s19);
        ipp=1/str2num(s19(4:ls19));
        
        ls22=length(s22);
        pw=str2num(s22(13:ls22))/1e9;
        
%        ls24=length(s24);
%        delay1=str2num(s24(16:ls24))/1e9;
        
        ls11=length(s11);
        nrg1=str2num(s11(6:ls11));
        
        ls13=length(s13);
        cint=str2num(s13(13:ls13));
        
        ls14=length(s14);
        r_interval=str2num(s14(9:ls14));
        
        ls2=length(s2);
        ns2=str2num(s2(9:ls2));
            
 %        ntrial=11600/cint/ipp/(nrg1)/(nch*2);
%        ntrial=24000;
%        if ns2/ipp/cint > 1024 & ns2/ipp/cint <2048
%        ntrial=fix(ns2/ipp/cint)+24;
        ntrial=fix(ns2/ipp/cint);
%        ntrial=512;
        
%        pause
%        n20=size(s23);
%        if n20(2)==33
%        antcode=[s23(32),s20(33)];
%        elseif n20(2)==32
%        antcode=[s23(32)];
%        end
   
       ls28=length(s28);
       min_range=str2num(s28(6:ls28))/1000;
        
        fprintf(fod,'date %s\r\n',year_date_time);
        fprintf(fod1,'%d ',year0);
        fprintf(fod1,'%d ',mon0);
        fprintf(fod1,'%d\r\n',day0);
        fprintf(fod,'time %s\r\n',num2str(time1));
        fprintf(fod,'beam_direction %s %s\r\n',s3_0(16:ls3_0),s3_1(ls3_1));
        fprintf(fod1,'%d ',h0);
        fprintf(fod,'hour %s\r\n',sh0);
        fprintf(fod1,'%d ',m0);
        fprintf(fod,'min  %s\r\n',sm0);
        fprintf(fod1,'%d',s0);
        fprintf(fod,'sec  %s\r\n',ss0);
        fprintf(fod,'nch  %d\r\n',nch);
        fprintf(fod,'prf  %s\r\n',s19(5:ls19));
        fprintf(fod,'ipp  %f\r\n',ipp);
        fprintf(fod,'pw   %f\r\n',pw );
        fprintf(fod,'min_range %d\r\n',min_range);
%        fprintf(fod,'delay1 %d\r\n',delay1);
        fprintf(fod,'nrg1 %d\r\n',nrg1);
        fprintf(fod,'cint %d\r\n',cint);
        fprintf(fod,'ntrial %d\r\n',ntrial);
        fprintf(fod,'ACQ_Time %d\r\n',ns2);
        %fprintf(fod,'antcode %s\r\n',antcode);
    %    fclose(fid);
        fclose(fod);
        fclose(fod1);
