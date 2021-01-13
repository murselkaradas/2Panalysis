## Copied from  Spontaneous behaviors drive multidimensional, brain-wide neural activity
*Carsen Stringer*1;2, Marius Pachitariu*1;3;4, Nicholas Steinmetz5, Charu Reddy5, Matteo Carandiniy5 and
Kenneth D. Harrisy3;4*

The matrices are too large to decompose in their raw form.
We therefore computed the SVD in two stages: 
first we  computed the SVD of short temporal segments, then we
concatenated the SVDs from different segments, scaled
by the singular values, and recomputed their SVD. Each
segment of frames is a matrix $G_i$. Since the number of
pixels is very large (> 1 million), we avoid computing the
SVD of this matrix directly, and instead compute the time
by time covariance matrix $G^T_iG_i$. We keep the top 200 eigenvectors $V_i$ of this matrix, which are also the top 200 right singular vectors of $G_i$. We then compute the spatial projections of these components 

$$U_i = G_iV_i$$ 

Notice that $U_i$ consists of the left singular vectors, scaled by the singular values. As such, the matrix  $U_i$ is a 200-dimensional summary of the segment $G_i$. We then concatenated $U_i$ for all segments of the movie, and re-computed the SVD:

$$[U \ S \ V^T ] = svd ([U_1 \ : \ : \  : \  U_n])$$
We kept the top **1,000** components of this $U$ matrix as the
spatial components of the face motion. We then projected
the raw movies onto these spatial components, to obtain
their temporal profiles:

$$W_{motion} = U^T G$$


### Compute fast SVD
```
function [U,S,V] = svdecon(X)

    [m,n] = size(X);

    if  m <= n
        C = X*X';
        [U,D] = eig(C);
        clear C;
    
        [d,ix] = sort(abs(diag(D)),'descend');
        U = U(:,ix);    
    
        if nargout > 2
            V = X'*U;
            s = sqrt(d);
            V = bsxfun(@(x,c)x./c, V, s');
            S = diag(s);
        end
    else
        C = X'*X; 
        [V,D] = eig(C);
        clear C;
    
        [d,ix] = sort(abs(diag(D)),'descend');
        V = V(:,ix);    
    
        U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
        %s = sqrt(sum(U.^2,1))';
        s = sqrt(d);
        U = bsxfun(@(x,c)x./c, U, s');
        S = diag(s);
end 
```


## split_SVD_2P for calcium imaging

The input parameter using the function ```splitSVD_2P.m```. 
| Parameter name | Description |
|----------------|-------------|
| ```frames``` | Nx $\times$ Ny $\times$ Nframes |
| ```num_block``` | number of block to divide Nframe into small set|
| ```num_sval```| number of singular top singular values kept |

Fsub : Nx $\times$ Ny $\times$ Nframe/num_block

For each block do following
```
[Ub,Sb,~] = svd(subF,'econ');
G{i}=Ub(:,1:num_svals)*Sb(1:num_svals,1:num_svals);
```
merge G and calculate SVD
```
G_all=cell2mat(G);
[U,svals,~] = svd(G_all,'econ');
SV=single(frames)'*U;
```

| Output parameter name   | Description |
|-------------------------|-------------|
| ```U``` | left singular vectors, Nx Ny $\times$ num_sval |
| ```SV``` | Projection of frame into left singular vectors U, Nframes $\times$ num_sval|
| ```sval```| top num_sval singular values |
