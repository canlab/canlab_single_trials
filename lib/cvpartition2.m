% same as cvpartition only it allows for respecting certain groupings (like
% a subject block) across partitions.
%
% Input ::
%
%   See 'help cvpartition' for most options
%
% Optional Input ::
%
%   'Stratify'  - Followed by vector with one element per observation
%                   containing a label specifying that observations
%                   grouping. For instance if you have 5 observations from
%                   10 subjects this vector might be something like,
%                   [1,1,1,1,1,2,2,2,2,2,3,3, ...,9,9,9,9,9,10,10,10,10,10]
%
%
% Written by Bogdan Petre, sometime 2019
classdef cvpartition2 < cvpartition
    properties (SetAccess = protected)
        sid;
    end
    methods
        % C = cvpartition2(GROUP, 'KFOLD',K,'Stratify',sid)
        function cv = cvpartition2(varargin)
            if ~isvector(varargin{5}) || length(varargin{1}) ~= length(varargin{5})
                error('sid must be a vector of sid''s of length(GROUP)');
            end
            
            sid = [];
            delete = [];
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'Stratify'
                            sid = varargin{5};
                            [~,b] = unique(sid);
                            varargin{1} = varargin{1}(b);
                            delete = i:i+1;
                    end
                end
            end
            varargin(delete) = [];
            
            cv@cvpartition(varargin{:});
            Impl = cvpartitionInMemoryImpl2(sid,varargin{:});
            Impl = Impl.updateParams();
            cv.Impl = Impl;
            
            cv.sid = sid;
        end % cvpartition constructor
        
        %{
        function testidx = test(cv,varargin)
            testidx = test(cv.Impl,varargin{:});
            uniq_sid = unique(cv.sid);
            test_sid = uniq_sid(testidx);
            testidx = ismember(cv.sid,test_sid);
        end
        
        function trainidx = training(cv,varargin)
            trainidx = training(cv.Impl,varargin{:});
            uniq_sid = unique(cv.sid);
            train_sid = uniq_sid(trainidx);
            trainidx = ismember(cv.sid,train_sid);
        end
        %}
    end    
end