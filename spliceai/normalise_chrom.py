
def normalise_chrom(source, target):
    ''' standardise chrom str to match a target style

    Args:
        source: chromsome string to standardise
        target: example of chromosome string style to match against

    Returns:
        chromosome string with 'chr' prefix added or removed as required
    '''
    has_prefix = lambda x: x.startswith('chr')
    
    if has_prefix(source) and not has_prefix(target):
        return source.strip('chr')
    elif not has_prefix(source) and has_prefix(target):
        return 'chr' + source
    return source
