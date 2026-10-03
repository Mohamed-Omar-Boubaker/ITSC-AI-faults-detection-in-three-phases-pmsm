classdef ChannelFeature < WindowFeature
    % Applies a FeatureGroup to one channel of the window ('IPos','INeg','Is',...)
    properties (SetAccess = immutable)
        Prefix
        Channel
        Group
    end
    methods
        function obj = ChannelFeature(prefix, channel, group)
            obj.Prefix = prefix; obj.Channel = channel; obj.Group = group;
        end
        function fs = extract(obj, window)
            fs = obj.Group.extract(window.(obj.Channel)).withPrefix(obj.Prefix);
        end
    end
end
