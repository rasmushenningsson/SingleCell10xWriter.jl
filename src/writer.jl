function string_dtype(len::Integer; fixed::Bool=true, utf8::Bool=false)
	dtype = datatype("")
	if fixed
		HDF5.API.h5t_set_size(dtype, max(len,1))
		HDF5.API.h5t_set_strpad(dtype, HDF5.API.H5T_STR_NULLPAD)
	else
		HDF5.API.h5t_set_size(dtype, HDF5.API.H5T_VARIABLE)
		HDF5.API.h5t_set_strpad(dtype, HDF5.API.H5T_STR_NULLTERM)
	end
	HDF5.API.h5t_set_cset(dtype, utf8 ? HDF5.API.H5T_CSET_UTF8 : HDF5.API.H5T_CSET_ASCII)
	dtype
end
string_dtype(str; kwargs...) = string_dtype(length(codeunits(str)); kwargs...)


function strings2fixed(v,len)
	mem = zeros(UInt8, length(v)*len)
	for (i,s) in enumerate(v)
		k = (i-1)*len+1
		c = codeunits(s)
		l = length(c)
		@assert l<=len
		mem[ k:k+l-1 ] .= codeunits(s)
	end
	mem
end


# --- Low level writing ---
function write_vec(g, name, v; chunk=80_000)
	chunk = min(length(v), chunk)
	g[name, chunk=(chunk,), shuffle=(), deflate=4] = v
	nothing
end

# use kwargs for chunk, shuffle and deflate
function write_strings(g, name, v; kwargs...)
	dtype = string_dtype(maximum(length∘codeunits,v))
	ds = create_dataset(g, name, dtype, size(v); kwargs...)
	mem = strings2fixed(v, sizeof(dtype))
	write_dataset(ds, dtype, mem)
end

function write_compressed_strings(g, name, v;
                                  chunk=10_000,
                                  filters=[Filters.Shuffle(), Filters.Deflate(4)],
                                  kwargs...)
	chunk = min(length(v), chunk)
	write_strings(g, name, v; chunk=(chunk,), filters, kwargs...)
end




function write_attr_string(g, name, value::String; kwargs...)
	dtype = string_dtype(value; kwargs...)
	attr = create_attribute(g, name, dtype, dataspace(value))
	write_attribute(attr, dtype, [value]) # if not vector it crashes...
end

function write_attr_string_vec(g, name, v::AbstractVector{<:AbstractString}; utf8::Bool=false)
	dtype = string_dtype(maximum(length∘codeunits,v); utf8)
	attr = create_attribute(g, name, dtype, dataspace(v))
	mem = strings2fixed(v, sizeof(dtype))
	# A hack, but I can't figure out how to write fixed string vectors for attributes in any other way
	mem2 = reinterpret(HDF5.FixedString{sizeof(dtype),0}, mem)
	write_attribute(attr, dtype, mem2)
end

write_attr(g, name, value) = attributes(g)[name] = value



# --- High level writing ---

function write_cellranger_attributes(h5; library_ids, original_gem_groups, version)
	write_attr_string(h5, "chemistry_description", "Single Cell 3' v3"; fixed=false)
	write_attr_string(h5, "filetype", "matrix"; utf8=true, fixed=false)
	write_attr_string_vec(h5, "library_ids", library_ids)
	write_attr(h5, "original_gem_groups", original_gem_groups)
	write_attr(h5, "version", version)
end


function write_cellranger_matrix(g, matrix)
	P,N = size(matrix)
	data = nonzeros(matrix)
	indices = rowvals(matrix) .- 1
	indptr = matrix.colptr .- 1

	# These conversions are rather silly, but we're mimicking an output file from cellranger
	# Might have been fixed in later cellranger versions though
	write_vec(g, "data", convert(Vector{Int32},data))
	write_vec(g, "indices", convert(Vector{Int64},indices))
	write_vec(g, "indptr", convert(Vector{Int64},indptr))
	write_vec(g, "shape", Int32[P,N])
end

function write_cellranger_barcodes(g, barcodes)
	write_compressed_strings(g, "barcodes", barcodes)
end

create_default(::Any, default::AbstractVector) = default
create_default(P, default) = fill(default,P)
get_or_create(df, key, default) = hasproperty(df,key) ? df[!,key] : create_default(size(df,1),default)

function write_cellranger_features(g, features; feature_type="Gene Expression", genome="hg19", read="", pattern="", sequence="")
	P = size(features,1)
	write_compressed_strings(g, "name", features.name)
	write_compressed_strings(g, "id", features.id)

	write_compressed_strings(g, "feature_type", get_or_create(features, "feature_type", feature_type))
	write_compressed_strings(g, "genome", get_or_create(features, "genome", genome))
	write_compressed_strings(g, "read", get_or_create(features, "read", read))
	write_compressed_strings(g, "pattern", get_or_create(features, "pattern", pattern))
	write_compressed_strings(g, "sequence", get_or_create(features, "sequence", sequence))
	
	write_strings(g, "_all_tag_keys", ["genome", "read", "pattern", "sequence"])
end

function write_cellranger_h5(h5::HDF5.File, counts::DataMatrix{<:SparseMatrixCSC};
                             library_ids,
                             original_gem_groups = [1],
                             version = 2,
                             kwargs...)
	write_cellranger_attributes(h5; library_ids, original_gem_groups, version)

	gmatrix = create_group(h5, "matrix")
	write_cellranger_matrix(gmatrix, counts.matrix)
	write_cellranger_barcodes(gmatrix, counts.obs.barcode)
	gfeatures = create_group(gmatrix, "features")
	write_cellranger_features(gfeatures, counts.var; kwargs...)
end

write_cellranger_h5(filename, counts; kwargs...) = h5open(filename, "w") do h5
	write_cellranger_h5(h5, counts; kwargs...)
end
