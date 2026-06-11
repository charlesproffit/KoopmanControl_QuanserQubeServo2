classdef Logger < handle
    %LOGGER Class to log data-driven synthesis procedure
    % 
    %   ------------------------------------------------------------------------
    %   Copyright 2023 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %
    
    properties
        fileID = 1;
        n_props
        syntax
        err_syntax
        divider
        title
        
        iter_syntax = "| % 4d |";
        iter_printed = 0;
        
        printf = @(fid, varargin) arrayfun(@(f) fprintf(f,  varargin{:}), fid, "UniformOutput", 0);
    end
    
    methods
        function obj = Logger(fileID, prop_name, prop_format)
            arguments
                fileID
            end
            arguments (Repeating)
                prop_name string
                prop_format string
            end

            if fileID == 0
                obj.printf = @(fid, varargin) [];
            end
            
            obj.fileID = fileID;
            obj.n_props = length(prop_name);
            obj.syntax = sprintf(" %s |", prop_format{:}) + "\n";
            
            tmp = @(x) regexp(x, "(?:%(?:\d\$)?[-+ 0#]?)(\d+)[^a-z$.]?", "tokens");
            field_widths = cellfun(@(y) cellfun(@double, tmp(y)), prop_format);
            
            obj.divider = "|------|" + sprintf("%s|", arrayfun( ...
                @(y) string(repelem('-', y+2)), ...
                field_widths)) + "\n";
            
            val = sum(field_widths) + length(field_widths)*3 - 3;
            obj.err_syntax = " %-" + string(val) + "s |\n";
            
            obj.title = "| Iter |" + ...
                sprintf(sprintf(" %%-%ds |", field_widths), prop_name{:}) + ...
                "\n";
        end
        
        function init(obj)
            %INIT Initial print for the logging table
            %   Prints the header and the horizontal divider denoting the start
            %   of the log.
            % 
            obj.printf(obj.fileID, obj.divider);
            obj.printf(obj.fileID, obj.title);
            obj.printf(obj.fileID, obj.divider);
        end
        
        function cleanup(obj)
            %CLEANUP Cleanup for the logging table
            %   Prints the last horizontal divider denoting the end of the log.
            % 
            obj.printf(obj.fileID, obj.divider);
        end
        
        function log_iter(obj, iter)
            %LOG_ITER Prints the current interation number
            % 
            if obj.iter_printed ~= 0
                obj.printf(obj.fileID, repmat('\b', 1, obj.iter_printed));
            end
            obj.printf(obj.fileID, obj.iter_syntax, iter);
            obj.iter_printed = 1;
        end
        
        function log(obj, varargin)
            %LOG Prints the information about the current iteration step
            % 
            if obj.iter_printed == 0
                obj.log_iter(0);
            end
            if length(varargin) ~= obj.n_props
                warning("Incorrect number of data");
                return;
            end
            obj.printf(obj.fileID, obj.syntax, varargin{:});  
            obj.iter_printed = 0;
        end
        
        function log_err(obj, varargin)
            %LOG_ERR Prints the error in the current iteration step
            % 
            if obj.iter_printed == 0
                obj.log_iter(0);
            end
            if length(varargin) ~= 1
                warning("Only single error message is allowed.");
                return;
            end
            obj.printf(obj.fileID, obj.err_syntax, varargin{:});
            obj.iter_printed = 0;
        end
    end
end

