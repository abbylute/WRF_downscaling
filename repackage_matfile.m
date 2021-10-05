function [status] = repackage_matfile(filename_in, filename_out, varname, NUM_X,NUM_Y,intFormat, SCALE, ADDOFFSET)
%
%[status] = repackage_matfile('/data/abby/ACLWDNB_CTRL_210m_chunk99.mat','/data/abby/ACLWDNB_CTRL_210m_chunk99_integer.mat', 'datdown', 130375, 28488,'int16', 100.0);
%filename_in = '/data/abby/ACLWDNB_CTRL_210m_chunk99.mat';
%filename_out = '/data/abby/ACLWDNB_CTRL_210m_chunk99_int16.mat';
%varname = 'datdown';
%NUM_X = 130375;
%NUM_Y = 28488;
%intFormat = 'int16'; 		%-this is the final integer format to save the output file
%SCALE=100.0; 			%-this is a scale factor to divide out of the data before saving as integers
			% for example: SCALE = 0.1 means 1 decimal digit of the data is preserved
%ADDOFFSET=0.0;			%-this is an add offset to subtract from the data before saving as integers


        %================================
	status = 0;

        %================================
        %    FILL VALUE
        %================================
        if(strcmp(intFormat,'int32'))
                maxvalue = 2^31;
                minvalue = -maxvalue;
                fillvalue = int32(minvalue); %min possible for int32=2^31
        elseif(strcmp(intFormat,'int16'))
                maxvalue = 2^15;
                minvalue = -maxvalue;
                fillvalue = int16(minvalue); %min possible for int16 =2^15
        elseif(strcmp(intFormat,'int8'))
                maxvalue = 2^17;
                minvalue = -maxvalue;
                fillvalue = int8(minvalue); %min possible for int16 =2^7
        elseif(strcmp(intFormat,'uint16'))
                maxvalue = 2^16;
                minvalue = -maxvalue;
                fillvalue = uint16(minvalue);
        else
                fillvalue =-9999;
        end

	%========================
	% SUBSETTING
	%========================
	x_split = 5; y_split=5;
	x_range = [1:ceil(NUM_X/x_split):NUM_X NUM_X];
	y_range = [1:ceil(NUM_Y/y_split):NUM_Y NUM_Y];
    % may need to move this below "abby's files" and read in the data dims
    % to set num_x and num_y instead of having them as arguments to the
    % function
    
	%========================
	% ABBY'S FILES
	%========================
    filehandle_in=matfile(filename_in); 

    
    %================================
    %   FIND MIN/MAX VALUE IN DATA
    %================================
	data_min = 10000000;	% where do these numbers come from and what is the point of them?
	data_max = -10000000;	
	start_y = 1;
	for j=1:length(y_range)-1
		end_y = y_range(j+1);
		range_y = [start_y:end_y];
		start_x = 1;
		for i=1:length(x_range)-1
			end_x = x_range(i+1);
			range_x = [start_x:end_x];
			data = filehandle_in.(varname)(range_x,range_y); 

			this_min= min(data(:));
			this_max= max(data(:));
			if(this_min<data_min)
				data_min = this_min;
            end
			if(this_max>data_max)
				data_max = this_max;
            end

			start_x = end_x + 1;
        end
		start_y = end_y + 1;
    end

    
    %================================
    %    CHECK MIN/MAX TO PREVENT TRUNCATION WITH INT FORMAT
    %================================
	if(data_min/SCALE<minvalue)
		disp('DATA_MIN < MINVALUE... You should change to a different intformat')
    end
	if(data_max/SCALE>maxvalue)
		disp('DATA_MAX >MAXVALUE... You should change to a different intformat')
    end
    
    
    %========================
    % NEW INTEGER FILE
    %========================
    filehandle_out = matfile(filename_out,'Writable',true);
    filehandle_out.(varname) = -fillvalue*ones(NUM_X,NUM_Y,intFormat);

    
        %================================
        %   FILL NEW MATFILE WITH INTEGERS INSTEAD OF SINGLES
        %================================
        start_y = 1;
        for j=1:length(y_range)-1
            end_y = y_range(j+1);
            range_y = [start_y:end_y];
            start_x = 1;
            for i=1:length(x_range)-1
                end_x = x_range(i+1);
                range_x = [start_x:end_x];

                data1 = filehandle_in.(varname)(range_x,range_y);
                data = round((data1-ADDOFFSET)/SCALE);
                data(isnan(data1)) = fillvalue;
                if(strcmp(intFormat,'int32'))
                    data_int = int32(data);
                elseif(strcmp(intFormat,'int16'))
                    data_int = int16(data);
                elseif(strcmp(intFormat,'int8'))
                    data_int = int8(data);
                elseif(strcmp(intFormat,'uint16'))
                    data_int = uint16(data);
                end

                filehandle_out.(varname)(range_x,range_y) = data_int;

                start_x = end_x + 1;
            end
            start_y = end_y + 1;
        end
	filehandle_out.SCALE = SCALE;
	filehandle_out.ADDOFFSET = ADDOFFSET;
	filehandle_out.directions=' VALUES = SCALE * INTEGERS_IN_FILE + ADDOFFSET';
        %================================
	success = 1;
end
