#####################################################################################################################
##
## The module is part of gkutils framework, for processing structures& sequences of biological macromolecules.
## It may not be copied without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2019)
##
#####################################################################################################################
from __future__ import print_function
import sys, os, os.path, types, gzip, json
import itertools
############################################################################################

def first(it):
    for x in it:
        return x
    raise ValueError('Empty list')
    
def shQuote(fn):
    r = ''
    for l in fn:
        if l in '()[]':
            r+=chr(92) +l
        else:
            r+=l
    return r
# Iterates through an iterable yielding chunked results
def iterchunks(lyst, length):
    r = []
    for v in lyst:
        r.append(v)
        if len(r)>=length:
            yield r
            r=[]
    if r:
        yield r

# Flatten the argument
def flatten(*it):
    for s in it:
        for sv in s:
            yield sv

# Dummy object corresponding to C struct construct
class StructObject(object):
    def __init__(self, **dictargs):    
        for k, v in list(dictargs.items()):
            object.__setattr__(self,k, v)
    def __getstate__(self): return self.__dict__
    def __setstate__(self, d): self.__dict__.update(d)     
    def __setattr__(self, k, v):
        object.__setattr__(self, k, v)
    def __str__(self):
        return "( " + ", ".join(["%s: %s" % (k, self.__dict__[k]) for k in sorted(self.__dict__)]) + " )";
    def __repr__(self):
        return self.__str__()
    
class UndefStructObject(StructObject):
    def __getattr__(self, k):
        return self.__dict__.get(k, None)
    def __getstate__(self): return self.__dict__
    def __setstate__(self, d): self.__dict__.update(d)

# Utility function to convert a dictionary (JSON/JSON-line) to a hierarchy of UndefStructObjects with (Struct with undefined default attribute value)
def parseJsonStruct(indata, sanitize=False):
    if sanitize:
        def _sanitize(k):
            r = ''
            for l in k:
                if l.isalpha():
                    r+=l
                else:
                    r+='_'
            return r.lower()
    else:
        def _sanitize(k):
            return k
    if isinstance(indata, dict):
        return UndefStructObject( **dict( (_sanitize(k), parseJsonStruct(v, sanitize=sanitize) ) for k,v in indata.items()) )
    elif isinstance(indata,list):
        return [ parseJsonStruct(v, sanitize=sanitize) for v in indata ]
    else:
        try:
            f = float(indata)
            try:
                i = int(indata)
            except:
                i = None
            if i==f:
                return i
            else:
                return f
        except:
            return indata

def condGzipOpener(fn, mode='r'):
    if fn.endswith('.gz'):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)
    
def handlify(mode='r', close=False, is_method=False, opener_f = open):
    ''' Decorator on a parse/read function.method, so that if a string is supplied it will attempt to open a file
    and use it instead '''
    # Intermediate function is needed so that decorator creation can be customised, by options above
    def intermediate(func):
        if is_method:
            def result(self, fh, *args, **kwargs):
                close_ = close
                r = None
                # basestring-->str in anticipation of python3
                if isinstance(fh, str):
                    fh = opener_f(fh, mode)
                    close_ = True
                try:                    
                    r = func(self, fh, *args, **kwargs)
                finally:
                    # Can't close if it's a generator, then the damn thing would loose the handle
                    if close_ and not isinstance(r, types.GeneratorType):
                        fh.close()
                return r
            return result
        else:                
            def result(fh, *args, **kwargs):
                close_ = close
                r = None
                # basestring-->str in anticipation of python3
                if isinstance(fh, str):
                    fh = opener_f(fh, mode)
                    close_ = True
                try:
                    r = func(fh, *args, **kwargs)
                finally:
                    if close_ and not isinstance(r, types.GeneratorType):
                        fh.close()
                return r
        return result
    return intermediate

from collections import defaultdict

#def _innerDefaultMaker(t):
#    return t

class __innerDefaultFactory(object):
    def __init__(self, t):
        self.t = t
    def __call__(self):
        return defaultdict(self.t)
    
def innerDefaultdict(t=None):
    #def __inner():
    #    return defaultdict(t)
    return __innerDefaultFactory(t)

from contextlib import contextmanager
@contextmanager
def changedDir(new_dir):
    try:
        old_dir = os.getcwd()
        os.chdir(new_dir)
        yield new_dir
    finally:
        os.chdir(old_dir)

import tempfile, shutil
@contextmanager
def tempDir(suffix='.tmpdir', prefix='tmp.', dir='.', change_dir=True):
    try:
        old_dir = os.getcwd()
        tmp_dir = tempfile.mkdtemp(suffix, prefix, dir)
        if change_dir:
            os.chdir(tmp_dir)
        yield tmp_dir
    finally:
        if change_dir:
            os.chdir(old_dir)
        shutil.rmtree(tmp_dir)

# Helper function to provide a temporary directory to a function
import tempfile, shutil, types
import contextlib
@contextlib.contextmanager
def tmpfile(give_name=False, ext='.tmp', prefix='tmp', dir=None):
    try:
        tmp_fh, tmp_fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
        if give_name:
            #os.close(tmp_fh)
            yield tmp_fh, tmp_fname
        else:
            yield tmp_fh
    finally:
        try:
            os.close(tmp_fh)
        except:
            pass
        try:
            os.unlink(tmp_fname)
        except:
            pass

@contextlib.contextmanager
def tmpdir(suffix='.tmp', prefix='tmp', dir='/tmp'):
    try:
        tmp_dir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
        yield tmp_dir
    finally:
        try:
            shutil.rmtree(tmp_dir)
        except:
            pass

def provideTempdir(suffix='.tmp', prefix='tmp', dir='/tmp', delete_dir=True, change_dir=False, is_method=False, is_generator=False):
    def intermediate(func):
        if not is_method:
            if not is_generator:
                def result(*args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        return func(dirname, *args, **kwargs)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
            else:
                def result(*args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        for r in func(dirname, *args, **kwargs):
                            yield r
                        if change_dir:
                            os.chdir(old_dir)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
        else:
            if not is_generator:
                def result(self, *args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        return func(self, dirname, *args, **kwargs)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
            else:
                def result(self, *args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        for r in func(self, dirname, *args, **kwargs):
                            yield r
                        if change_dir:
                            os.chdir(old_dir)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
        return result
    return intermediate

def provideTempfile(give_name=False, ext='.tmp', prefix='tmp', dir=None, is_method=False):
    def intermediate(func):
        if not is_method:
            def result(*args, **kwargs):
                fh, fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
                try:
                    if give_name:
                        return func(fh, fname, *args, **kwargs)
                    else:
                        return func(fh, *args, **kwargs)
                finally:
                    try:
                        os.close(fh)
                    except:            
                        pass
                    try:
                        os.unlink(fname)
                    except:
                        pass
        else:
            def result(self, *args, **kwargs):
                fh, fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
                try:
                    if give_name:
                        return func(self, fh, fname, *args, **kwargs)
                    else:
                        return func(self, fh, *args, **kwargs)
                finally:
                    try:
                        os.close(fh)
                    except:            
                        pass
                    try:
                        os.unlink(fname)
                    except:
                        pass
        return result
    return intermediate
         
#######################################################################################################################
## Parsing/writing basic forms of list, set and dictionary (1/line, dicts in tab/comma/whatever-separated files)
#######################################################################################################################
@handlify(is_method=False)
def parseSimpleList(fh, type_=str):
    result = list()
    for line in fh:
        line = line.strip()
        if line:
            result.append( type_(line) )
    return result

@handlify(mode='w', is_method=False)
def writeSimpleList(fh, lyst):
    for v in lyst:
        print("%s" % (v), file=fh)

writeSimpleSet = writeSimpleList

@handlify(is_method=False)
def parseSimpleSet(fh, type_=str):
    result = set()
    for line in fh:
        line = line.strip()
        if line:
            result.add( type_(line) )
    return result

@handlify(is_method=False)
def parseSimpleDict(fh, separator='\t', key_f=None, value_f=None, none_val='None'):
    ''' Parses a dictionary in form of multiple character-separated values (one line, one key-value pair)'''    
    result = {}
    for line in fh:
        line = line.strip()
        if not line:
            continue
        fs = line.split(separator, 1)
        try:
            k, v = fs[0], fs[1]
        except:
            raise ValueError(line)
        if key_f is not None:
            k = key_f(k)
        if value_f is not None:
            v = value_f(v)
        if none_val and none_val==v:
            v = None
        result[k] = v
    return result

def reverseDict(d):
    ''' Reverses a dictionary '''
    return dict( (v,k) for k,v in d.items() )


@handlify(mode='w', is_method=False)
def writeSimpleDict(fh, dikt, separator='\t', sort=False, klower=False, **kwds):
    ''' Writes a dictionary in form of multiple character-separated values (one line, one key-value pair)'''
    #return None
    #print >> sys.stderr, "GUMBAS"
    if sort:
        for k,v in sorted(iter(dikt.items()), key=lambda v:v[0].lower() if klower else v[0]):
            #print k,v
            print("%s%s%s" % (k,separator,v), file=fh)
    else:
        for k,v in dikt.items():
            print("%s%s%s" % (k,separator,v), file=fh)
    if 'eof' in kwds:
        print(kwds['eof'], file=fh) 


def iterateDir(in_dir, prefix=None, suffix=None):
    ''' Process a directory into a set of filenames (optionally by suffix or prefix) '''
    for fname in os.listdir(in_dir):
        if (prefix is None or fname.startswith(prefix)) and (suffix is None or fname.endswith(suffix)):
            full_fname = os.path.abspath(os.path.join(in_dir, fname))
            yield full_fname

@contextmanager
def condMakeDir(out_dir):
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    with changedDir(out_dir):
        yield out_dir
        
def fname2pdb(fname):
    return os.path.basename(fname).split('.')[0].lower()

def stringize(it):
    ''' Convert an iterables\' elements to strings. '''
    for i in it:
        yield str(i)

@handlify(mode='r', is_method=False)
def endsWithDone(fh):
    done = False
    for line in fh:
        if line.endswith('#DONE'):
            done=True
            break
    return done

def existsDone(fname):
    ''' Does a .done marker file exists '''
    return os.path.exists(fname+'.done')

def markAsDone(fname):
    '''  Creates a .done marker file '''
    open(fname+'.done', 'w').close()

def hasContent(fname):
    return (os.path.exists(fname) and os.path.isfile(fname) and not os.path.getsize(fname)==0)

def parseResidueId(s):
    ''' Parse residue identifier into (number, insertion code) tuple.'''
    s = s.strip()
    if s[-1] in '0123456789':
        resno, icode = int(s), ' '
    else:
        resno, icode = int(s[:-1]), s[-1]
    return resno,icode

#######################################################################################################################
## Parsing/writing basic forms of list, set and dictionary (1/line, dicts in tab/comma/whatever-separated files)
#######################################################################################################################
@handlify(is_method=False)
def parseSimpleList(fh, type_=str):
    result = list()
    for line in fh:
        line = line.strip()
        if line:
            result.append( type_(line) )
    return result

@handlify(mode='w', is_method=False)
def writeSimpleList(fh, lyst):
    for v in lyst:
        print("%s" % (v), file=fh)

#@handlify(is_method=False)
#def parseSimpleSet(fh):
#    result = set()
#    for line in fh:
#        line = line.strip()
#        if line:
#            result.add(line)
#    return result

@handlify(is_method=False)
def parseSimpleDictX(fh, separator='\t', key_f=None, value_f=None):
    ''' Parses a dictionary in form of multiple character-separated values (one line, one key-value pair)'''    
    result = {}
    for line in fh:
        line = line.strip()
        fs = line.split(separator, 1)
        k, v = fs[0], fs[1]
        if key_f is not None:
            k = key_f(k)
        if value_f is not None:
            v = value_f(v)
        result[k] = v
    return result

#@handlify(mode='w', is_method=False)
#def writeSimpleDictX(fh, dikt, separator='\t', **kwds):
#    ''' Writes a dictionary in form of multiple character-separated values (one line, one key-value pair)'''
#    for k,v in dikt.iteritems():
#        print >> fh, "%s%s%s" % (k,separator,v)
#    if kwds.has_key('eof'):
#        print >> fh, kwds['eof'] 

@handlify(is_method=False)
def parseSetDict(fh, separator='\t', key_f=None, secsep=',', value_f=None):
    ''' Parses a dictionary in form of multiple character-separated values (one line, one key-value pair)'''    
    result = {}
    for line in fh:
        line = line.strip()
        fs = line.split(separator, 1)
        k, v = fs[0], fs[1]
        vs = v.split(secsep)
        if key_f is not None:
            k = key_f(k)
        if value_f is not None:
            vs = [ value_f(v_) for v_ in vs ]
        result[k] = set(vs)
    return result

#######################################################################################################################
## IniOptionParser - extended ConfigParser defaulting to .ini file for values
#######################################################################################################################
try:
    import configparser, optparse
except:
    import configparser, optparse

class IniOptionParser(optparse.OptionParser, object):
    ''' Subclassed OptionParser from optparse, which can take arguments from the parsed .INI config file. '''
    def __init__(self, inifile, *args, **kwargs):
        optparse.OptionParser.__init__(self, *args, **kwargs)
        self.add_option("--inifile", dest="inifile",
                        default=inifile,
                        help="specify a file containing further specified options (if invalid print a warning and proceed) [%s]" % (inifile)
                        )
    def _dictify(self, ini_opts, types={}):
        result = {}
        for section in ini_opts.sections():
            if section=='__MAIN__':
                inner = result
            else:
                inner = result[section] = {}
            for option, value in ini_opts.items(section):
                if section in types and option in types[section]:
                    inner[option] = types[section][option](value)
                else:
                    inner[option] = value
        return result
    def parse_args(self, types={}):
        (options, args) = super(IniOptionParser, self).parse_args()
        if os.path.exists(options.inifile):
            parser = configparser.SafeConfigParser()
            parser.read([options.inifile])
            ini_options = self._dictify(parser, types)
            for key, value in ini_options.items():
                self.defaults[key] = value
            (options, args) = super(IniOptionParser, self).parse_args()
        return (options, args)

def makeListTokenizer(separator=",", create_f=str):
    def inner(data):
        return [create_f(x) for x in data.split(separator)]
    return inner

def separatingCallback(option, opt_str, value, parser, create_f=str, separator=","):
    if value:
        setattr(parser.values, option.dest, makeListTokenizer(separator, create_f)(value))

EXTENSIONS = {
    'vdw': '.residue.cs',
    'loops': '.loops',
    'locks': '.locks',
    'domains': '.domains'
    }
 

def quotev(v, qchar="'"):
    if v is None:
        return 'NULL'
    elif (isinstance(v, str) or isinstance(v, str)) and not (v.startswith(qchar) and v.endswith(qchar)):
        if isinstance(v, str):
            v = str(v)
        if "'" in v:
            v = v.replace("'", "\\'")
            #return '"'+v+'"'
        return str("'")+v+str("'")
    else:
        return v

def quoteInStmt(stmt, **vals):
    for k,v in vals.items():
        #try:
        vals[k] = quotev(v)
        #except:
        #    print "key", k
        #    print "val", v
        #    raise
    return stmt.format(**vals)

def formatForCopy(sep, *vals):
    #for k,v in vals.iteritems():
    #    vals[k] = unicode(v)
    return str(sep).join([str(v) for v in vals])+"\n"

#######################################################################################################################################
#
# Writing/reading matrices
#
#######################################################################################################################################

@handlify(mode='w')
def writeMatrixD(fh, mdikt, ordering):
    sep = '\t'
    for sid1 in ordering:
        print(sep.join( str(mdikt[sid1][sid2]) for sid2 in ordering), file=fh)
            
@handlify(mode='w')
def writeAdjListD(fh, mdikt, ordering):
    sep = '\t'
    seen = set()
    for sid1 in ordering:
        for sid2 in ordering:
            if not sid2 in seen:
                d = mdikt[sid1][sid2]
                print(sep.join( str(v) for v in [ sid1, sid2, d] ), file=fh)
            seen.add(sid2)
        
# Fixed so it defaults to text instead of binary            
def conditionalGzipOpen(fname, *params, **kwparams):
    #print("KUrwa", params, kwparams)
    if fname.endswith('.gz'):
        if 'mode' in kwparams and kwparams['mode'] in ('r','w'):
            kwparams['mode']+='t'
        elif params[0] in ('r','w'):
            params = list(params)
            params[0]=params[0]+'t'
        return gzip.open(fname, *params, **kwparams)
    else:
        return open(fname, *params, **kwparams)


import hashlib
def seqHash(txt):
    h = hashlib.md5()
    h.update(txt)
    return h.hexdigest()

####################################################
# Extract sid from eXtended Identifier (xid)
def sidFromXid(xid):
    if type(xid) is int:
        return xid
    if xid.endswith('_'):
        xid = xid.rstrip('_')
    try:
        sid = int(xid.split('__')[0])
        return sid
    except:
        return None
        
def seqidFromXid(xid):
    fields = xid.split('__',1)
    dbname = fields[0]
    seqid = fields[1]
    if seqid.startswith('_'):
        dbname=dbname+'_'
        seqid = seqid[1:]
    return StructObject(dbname=dbname, seqid=seqid)
#####################################################
# Parse a list with header (TSV format)
@handlify(mode='r', opener_f=conditionalGzipOpen)
def parseTsvList(fh, fnames=None, sep='\t', none_val='None', none_val_inner=None, include_line_attr=False, strip_fname_spaces=True, comment_chars='##'):
    def _valueMake(v):
        if none_val is not None and v==none_val:
            return none_val_inner
        elif v.endswith('"') and v.startswith('"'):
            return v.strip('"')
        else:
            try:
                return int(v)
            except:
                try:
                    return float(v)
                except:
                    return v
    for lineno, line in enumerate(fh):
        fields = line.strip().split(sep)
        if line.startswith(comment_chars):
            continue
        elif not fnames:
            if line.startswith('#'):
                fields = line.strip('#').strip().split(sep)
            fnames = fields
            if strip_fname_spaces:
                fnames = [ f.strip().replace(' ', '_').replace('.','_').replace('-','__').replace('(','_').replace(')','_').replace('?','_').replace('"','').replace("'",'') for f in fnames ]                
        elif line.startswith('#'):
            continue
        else:
            fields = [ _valueMake(v) for v in fields ]
            if len(fields)<len(fnames):
                delta = len(fnames) - len(fields)
                fields.extend( [ none_val_inner ] * delta )
            x = StructObject(**dict(list(zip(fnames, fields))))
            if include_line_attr:
                x.line = line.strip()
            if not hasattr(x, 'fnames'):
                x.fnames = fnames
            yield x

@handlify(mode='wt')
def writeTsvList(fh, recs, fnames=None, sep='\t', none_val='None', to_val_func=None, skip_fields=['fnames', 'line'], use_utf8=False):
    def _toValFunc(v):
        if v is None and none_val:
            return none_val
        if use_utf8:
            return str(v)
        else:
            return str(v)
    if to_val_func:
        _toValFunc = to_val_func
    for i, rec in enumerate(recs):
        if i==0:
            if not fnames:
                if 'fnames' in skip_fields and hasattr(rec, 'fnames'):
                    fnames = rec.fnames
                else:
                    fnames = sorted(k for k in list(rec.__dict__.keys()) if k not in skip_fields)
            print(sep.join(fnames), file=fh)
        # Not needed anymore - 17/05/2020 - different handling in Py3
        #if 0 and use_utf8:
        #    print("\t".join( [ _toValFunc(getattr(rec,fname, None)) for fname in fnames] ).encode('utf8'), file=fh)
        #else:
        print("\t".join( [ _toValFunc(getattr(rec,fname, None)) for fname in fnames] ), file=fh)
            
def toValue(v):
    try:
        return eval(v)
    except:
        return str(v)

def parseLabelDiktTxt(txt, label_sep=';', indiv_sep=',', label_div='=', reverse=True):
    txt = txt.strip()
    labtxts = txt.split(label_sep)
    r = defaultdict(list)
    for labtxt in labtxts:
        lab, vs = labtxt.split(label_div)
        vs = vs.split(indiv_sep)
        r[lab].extend(list(v.strip() for v in vs))
    if reverse:
        r_new = {}
        for lab in r:
            for v in r[lab]:
                r_new[v]=lab
        return r_new
    else:
        return r

# StructObject-aware JSON parsing
class StructJsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, StructObject):
            return obj.__dict__
        elif isinstance(obj, set) or isinstance(obj, frozenset):
            try:
                return sorted(obj)
            except:
                return obj
        return json.JSONEncoder.default(self, obj)
    
class StructJsonDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_pairs_hook=self.object_pairs_hook, *args, **kwargs)
    def object_pairs_hook(self, pairs):
        def _procPair(k,v):
            if k.endswith('_dikt') or k.endswith('_map'):
                return (k, v.__dict__)
            elif k.endswith('_set'):
                return (k, set(*v))
            return (k, v)
        return UndefStructObject(**dict( _procPair(k,v) for k,v in pairs ))
    
@handlify(mode='r')
def readJsonStruct(fh):
    return json.load(fh, cls=StructJsonDecoder)

@handlify(mode='w')
def writeJsonStruct(fh, v, indent=4):
    json.dump( v, fh, cls=StructJsonEncoder, indent=indent)
    fh.write('\n')
    

# CONSTANTS
try:
    import psutil
    NCPU = psutil.cpu_count(logical=False)
except:
    import multiprocessing
    NCPU = multiprocessing.cpu_count()
#NCPU = 4


