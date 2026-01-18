function err = fitHist(fname, param1, err1, verbose, time)
% function err = fitHist(fname, param, err, verbose, time)
% 
% remember the relationship between parameters and errors and return the error
% when the same parameters are used.
%
% usage:
%  %when remembering error
%  fitHist(fname, param, err);
%  %when recalling error
%  err = fitHist(fname, param);
%  %err is empty when no history.
%
%  %to reset (delete history)
%  fitHist(fname, 'reset')
%
%  %to set seed
%  fitHist(fname, 'seed');
%  or
%  fitHist(fname, 'seed', seed);
%  % if fname has seed, it will reset rand using that seed. If not create a new one.

if isempty(fname)
    err = [];
    return;
end

if ~exist('verbose', 'var')
    verbose = false;
end

if ~exist('time', 'var')
    time = [];
end

if ischar(param1) && strcmp(param1, 'reset')
    if exist(fname, 'file')
        delete(fname);
    end
    fprintf('[fitHist] Reset file: %s.\n', fname);
    err = [];
    return;
elseif ischar(param1) && strcmp(param1, 'seed')
    if exist(fname, 'file')
        try
            S = load(fname);
            if isfield(S, 'seed')
                seed = S.seed;
                if exist('err1', 'var') && ~isempty(err1) && err1 ~= seed
                    error('Specified seed does not match with the seed in file');
                end
            else
                error('seed could not be found in %s', fname);
            end
        catch
            warning('File %s seems broken. Starting from scratch.', fname);
            if exist('err1', 'var') && ~isempty(err1)
                seed = err1;
            else
                seed = round(sum(clock*1000));
            end
            save(fname, 'seed');
        end
    else
        if exist('err1', 'var') && ~isempty(err1)
            seed = err1;
        else
            seed = round(sum(clock*1000));
        end
        save(fname, 'seed');
    end
    fprintf('[fitHist] Set seed = %d.\n', seed);
    try
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
    catch
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
    end
    err = [];
    return;
elseif ischar(param1) && strcmp(param1, 'comment')
    if ~exist('err1', 'var') || isempty(err1)
        error('comment require the third input variable');
    end
    
    comment = err1; 
    if exist(fname, 'file')
        save(fname, 'comment', '-APPEND');
    else
        save(fname, 'comment');
    end
    err = [];
    return;
end




if exist('err1', 'var') && ~isempty(err1) % remember mode
    if ~exist(fname, 'file')
        param = param1(:)';
        err = err1;
        fit_time = time;
        save(fname, 'param', 'err', 'fit_time');
        if verbose
            fprintf('[fitHist] Saved parameters in %s (N = %d; err = %1.3f).\n', fname, 1, err1);
        end
    else
        S = load(fname);
        if ~isfield(S, 'param')
            param = param1(:)';
            err = err1;
            fit_time = time;
            save(fname, 'param', 'err', 'fit_time', '-APPEND');
            if verbose
                fprintf('[fitHist] Saved parameters in %s (N = %d; err = %1.3f).\n', fname, 1, err1);
            end
        else
            idx = ismember(S.param, param1(:)', 'rows');
            if any(idx)
                err2 = S.err(idx);
                % if abs(err1-err2)<1e-10 || (isnan(err1) && isnan(err2))
                if abs(err1-err2)<1e-6 || (isnan(err1) && isnan(err2)) % debug, Jiahao
                    fprintf('[fitHist] Parameters already exist.\n');
                    return; % already exists
                else
                    error('Inconsistency in errors: registerd %f, new %f', err2, err1);
                end
            end
            param = [S.param; param1(:)'];
            err = [S.err; err1];
            fit_time = [S.fit_time; time];
            save(fname, 'param', 'err', 'fit_time', '-APPEND');
            if verbose
                fprintf('[fitHist] Saved parameters in %s (N = %d; err = %1.3f).\n', fname, length(err), err1);
            end
        end
    end
    err = err(end);
else  % recall mode
    if ~exist(fname, 'file')
        err = [];
    else
        S = load(fname);
        if ~isfield(S, 'param')
            err = [];
        else
            idx = ismember(S.param, param1(:)', 'rows');
            if sum(idx)>1
                error('duplicated parameters');
            end
            err = S.err(idx);
        end
    end
    if verbose
        if isempty(err)
            fprintf('[fitHist] No parameters. Returned empty.\n');
        else
            [~, loc] = ismember(err, S.err);
            fprintf('[fitHist] Returned error %1.3f from %s (id = %d).\n', err, fname, loc); % loc ..where is err in S.err
        end
    end
end
























