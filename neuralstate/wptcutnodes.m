function [wpTree] = wptcutnodes(wpTree, nodes)
%WPTCUTNODES Cut nodes from the giving wavelet pocket tree.
%   Use as:
%       wpTree = wptcutnodes(wpTree, nodes);
%
%   Input:
%       wpTree  - the original wavelet pocket tree that need to be cutted
%       nodes   - {[x1,y1], [x2,y2], ...}, nodes to be cutted
%
%   Output:
%       wpTree  - cutted wavelet poecket tree

% Convert node index
nCutNodes = length(nodes);
cutNodesIdx = zeros(1, nCutNodes);
for iNode = 1:nCutNodes
    cutNodesIdx(iNode) = nodeidx(nodes{iNode});
end


% Cut nodes
nNodes = length(cutNodesIdx);
cutNodesIdx = sort(cutNodesIdx, 'descend');
for iNode = 1:nNodes
    cutNode = cutNodesIdx(iNode);
    
    cutLevl = nodeidx(cutNode);
    cutLevl = cutLevl(1);
    endNodes = nodeidx([cutLevl,0]):nodeidx([cutLevl,2^cutLevl]);
    wpTree = wpjoin(wpTree, endNodes);
        
    cfsCut = zeros(read(wpTree, 'sizes', cutNode));

    if mod(cutNode,2)==0
        cutNodePair = cutNode - 1;
    else
        cutNodePair = cutNode + 1;
    end
    cfsPair = wpcoef(wpTree, cutNodePair);

    wpTree = write(wpTree, 'cfs', cutNode, cfsCut, 'cfs', cutNodePair, cfsPair);

%     leftNode = nodeidx(min(cutNode, cutNodePair));
%     superNode = nodeidx([leftNode(1)-1, leftNode(2)/2]); % Super node - (level-1) and (index/2)
%     wpTree = wpjoin(wpTree, superNode);
end

end

