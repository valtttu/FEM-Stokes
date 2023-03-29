% Make a mesh composed of rectangular blocks based on a text file

% Parameters:
% path      Location for the domain defining text file
% N         Number of mesh refinements to be made, 0 => just two triangles
%           per block
% Returns:
% mesh      A mesh object that has the fields given by inittri + bdof & bvals
%           fields that allows the access to the boundary elements and
%           setting their values

function mesh = build_mesh(path, N)
    
    % Define keywords
    kws = {'points', 'blocks', 'bc'}; % the keywords to look from the file
    in.points = []; % a 2xNp matrix of the mesh points
    in.blocks = []; % a 4xNb matrix of the blocks
    in.bc = []; % a 4xNv matrix that contains the points defining the outer edge and ux and uy for that edge

    % Read the text file
    fid = fopen(path);
    while(~feof(fid))
        tline = fgetl(fid);
        if(tline == -1)
            break;
        end

        if(~contains(tline, '#'))
            % Found a keyword -> fill the corresponding field
            for i = 1:length(kws)
                if(strcmp(kws{i}, replace(tline, ' ', '')))
                    tline = fgetl(fid); % skip the bracket
                    tline = fgetl(fid);
                    j = 1;
                    while(~contains(tline, '}'))
                        in.(kws{i})(:,j) = str2num(tline)';
                        j = j + 1;
                        tline = fgetl(fid);
                    end
    
                end
            end
        end
    end

    fclose(fid);

    % Define the triangulation
    p = in.points;
    t = [];
    for i = 1:size(in.blocks, 2)
        t(:,2*i-1) = in.blocks([1,2,4],i) + 1; % +1 since the points are indexed from 0...
        t(:,2*i) = in.blocks([2,3,4],i) + 1;
    end

    % Make the mesh
    mesh = inittri(p, t);

    % Refine N-times
    for i=1:N
      mesh = refine_tri(mesh);
    end

    % Add the boundary conditions
    out_inds = find(mesh.e2t(2,:) == 0);
    out_edges = mesh.edges(:,out_inds);
    out_ts = unique(mesh.e2t(1,out_inds(:)));
    bdof = [];
    bvals = [];
    for i = 1:length(out_ts)
        % Go over the edges of triangle i
        ti = out_ts(i);
        for j = 1:3
            if(sum(mesh.t2e(j, ti) == out_inds) == 1)
                eind = mesh.t2e(j, ti) == out_inds;
                ei = j;
          
                % Append the set of all outer edge indices
                bdof(end+1) = mesh.edof(ti,ei,1);
                bdof(end+1) = mesh.edof(ti,ei,2);
        
                % Set the bounary values
                x = mean(mesh.p(1,out_edges(:,eind)));
                y = mean(mesh.p(2,out_edges(:,eind)));
                xl = max([p(1, in.bc(1,:)+1);p(1, in.bc(2,:)+1)]);
                xs = min([p(1, in.bc(1,:)+1);p(1, in.bc(2,:)+1)]);
                i1 = x - xs <= xl - xs & x - xs >= 0;
                yl = max([p(2, in.bc(1,:)+1);p(2, in.bc(2,:)+1)]);
                ys = min([p(2, in.bc(1,:)+1);p(2, in.bc(2,:)+1)]);
                i2 = y - ys <= yl - ys & y - ys >= 0; 
                ibc = i1 & i2;
                bvals(end+1) = in.bc(3,ibc);
                bvals(end+1) = in.bc(4,ibc);
            end
        end

        

    end

    mesh.bdof = bdof;
    mesh.bvals = bvals;
end