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
    end
    methods (Access = ?cvpartition2)        
        function obj = updateParams(obj)
            obj.N = length(obj.sid);
            for i = 1:obj.NumTestSets
                obj.TrainSize(i) = sum(ismember(obj.sid,find(obj.training(i))));
                obj.TestSize(i) = sum(ismember(obj.sid,find(obj.training(i))));
            end
        end
    end
end