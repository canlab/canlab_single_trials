% this method corrects the counters of cvpartitions, which under the
% cvpartition2 scheme don't count obsevations, and instead count blocks.
classdef cvpartitionInMemoryImpl2 < internal.stats.cvpartitionInMemoryImpl
    properties (SetAccess = private)
        sid = []
    end
    methods
        function obj = cvpartitionInMemoryImpl2(sid, varargin)
            obj = obj@internal.stats.cvpartitionInMemoryImpl(varargin{:});
            [~,~,obj.sid] = unique(sid,'stable');
        end
        
        function obj = repartition(obj, varargin)
            [~,b] = unique(obj.sid);
            obj.Group = obj.Group(b);
            obj.indices = obj.indices(b);
            obj.N = length(b);
            
            obj = repartition@internal.stats.cvpartitionInMemoryImpl(obj, varargin{:});
            
            obj = obj.updateParams();
        end
    end
    methods (Access = ?cvpartition2)        
        function obj = updateParams(obj)
            obj.N = length(obj.sid);
            for i = 1:obj.NumTestSets
                obj.TrainSize(i) = sum(ismember(obj.sid,find(obj.training(i))));
                obj.TestSize(i) = sum(ismember(obj.sid,find(obj.test(i))));
            end
            
            [newIndices, newGroup] = deal(zeros(obj.N,1));
            uniq_sid = unique(obj.sid);
            for i = 1:length(uniq_sid)
                this_sid = uniq_sid(i);
                sid_idx = find(this_sid == obj.sid);
                newGroup(sid_idx) = obj.Group(i);
                newIndices(sid_idx) = obj.indices(i);
            end
            obj.indices = newIndices;
            obj.Group = newGroup;
            
            if ~isempty(obj.holdoutT)
                warning(['Warning obj.holdoutT is not empty. Using this ',...
                    'function in this manner has not been tested. Please ',...
                    'check cvpartitionMemoryImpl2 and see how holdoutT is ',...
                    'being handled, compare with internal.stats.cvpartitionInMemoryImpl ',...
                    'and fixing it if necessary.']);
            end
        end
    end
end