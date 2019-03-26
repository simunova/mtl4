import getopt, sys, os, os.path, re, string

######################
# Helper functions
######################

# platform-dependent preprocessor definitions
def macro_flag(env, flag):
    if sys.platform == 'win32' and conf.env['CC'] == 'cl':
        return '/D' + flag
    else:
        return '-D' + flag    

# add debugging flags to environment
def add_debug(env):
    if sys.platform == 'win32' and conf.env['CC'] == 'cl':
        env.Append(CCFLAGS = '/Z7')
        env.Append(CCFLAGS = '/MDd')
    else:
        env.Append(CCFLAGS = '-g')
    #print macro_flag(env, 'MTL_ASSERT_FOR_THROW')
    env.Append(CCFLAGS = macro_flag(env, 'MTL_ASSERT_FOR_THROW'))
   

def check_g5():
    print "platform ", sys.platform
    if sys.platform == "darwin" :

	system("/usr/bin/system_profiler SPHardwareDataType")
	profile = os.popen("system_profiler SPHardwareDataType").read()
	print profile
	if profile.find("Mac G5") != -1:
	    print "is G5"
	    return 1
    print "isn't G5"
    return None


def add_platform_flags(env):
    if env['CC'] == 'cl':
        # Turn off some warnings: 1. for problems with CRTP 2. with VC's own headers 
        #                         3. alleged ambiguos assignment
        # I hope we find a cleaner way to get rid of this warning
        wd = ' /wd4355 /wd4996 /wd4522'
        env.Append(CCFLAGS = '/EHa /D_CONSOLE /D_CRT_SECURE_NO_DEPRECATE /DNOMINMAX' +
			     ' /D_SECURE_SCL=0' + wd)

        # The following two if shouldn't be necessary
        # Note that INCLUDE must be split and LIB must not :-P
        if os.environ.has_key('INCLUDE') :
            for p in os.environ['INCLUDE'].split(';'):
                env.Append(CPPPATH = p)
        if os.environ.has_key('LIB') :
            env.Append(LIBPATH = os.environ['LIB'])


# add optimization flags to environment
# normal optimization for opt=1, for opt=2 high optimization
# at the moment additional flags are added in this function
# -> this needs to be cleaned up!
def add_opt(env, opt):
    if int(opt):
        env.Append(CCFLAGS = macro_flag(env, 'NDEBUG'))
    if sys.platform == 'win32' and conf.env['CC'] == 'cl':
        if int(opt) == 1:
            env.Append(CCFLAGS = '/O2 /Ot /Ob2')
        elif int(opt) == 2:
            env.Append(CCFLAGS = '/Ox')
            env.Append(CCFLAGS = '/MD')
    else:
        if int(opt) == 2:
    	    env.Append(CCFLAGS = '-O3')
    	    if conf.env['CC'] == 'gcc':
    	        env.Append(CCFLAGS = '-ffast-math')
    	        #check_g5()
            # env.Append(CCFLAGS = '-pg') # for profiling
    	    # env.Append(LINKFLAGS = '-pg') # for profiling
        elif int(opt) == 1:
    	    env.Append(CCFLAGS = '-O2')
            # env.Append(CCFLAGS = '-pg') # for profiling
    	    # env.Append(LINKFLAGS = '-pg') # for profiling
        
    if env['CC'] != 'cl':
        env.Append(CXXFLAGS = conf.env['CCFLAGS'])
	

def check_no_long_double(conf):
    cc = conf.env['CC']
    if cc == 'gcc':
        # conf.env.Append(CXXFLAGS = '-Werror -Wreorder -pedantic -Wpointer-arith -Wcast-align -Wcast-qual -Wwrite-strings')
	output = os.popen(cc + " --version").read()
	st = output.split('\n')[0]
	tmatch = re.search(".* \(GCC\) (\d\.\d)", st)
	if tmatch and string.atof(tmatch.group(1)) == 4.0:
	    conf.env.Append(CCFLAGS = '-Wno-long-double')
	    # conf.env.Append(CXXFLAGS = '-Wno-long-double')


def add_user_opt(env, add_optflag):
     env.Append(CCFLAGS = add_optflag)
     env.Append(CXXFLAGS = add_optflag)

###############################
# check for installed blas ...
# search different libraries
###############################
def detected_blas(env): # blas library found
   env.Append(CPPDEFINES = 'MTL_HAS_BLAS')

def check_for_blas(env):
 # symbolname
 symname = 'dgemm_'
 found_blas = 0

 # get command line blas path
 blas_path = ARGUMENTS.get('blas_path', '')
 if blas_path:
   print 'adding ' + blas_path + ' to LIBPATH.'
   env.Append(LIBPATH = [ blas_path, blas_path + '/lib/' ])
 
 # extra linker flags for blas
 blas_ldflags = ARGUMENTS.get('blas_ldflags', '')
 if blas_ldflags:
   env.Append(_LIBFLAGS=blas_ldflags.split())

 # get command line blas lib
 blas_lib = ARGUMENTS.get('blas_lib', '')
 if blas_lib:
   myenv = env.Copy()
   conf = Configure(myenv)
   # check supplied lib
   print "Checking for lib " + blas_lib + "..."
   if conf.CheckLib(blas_lib, symname):
     myenv = conf.Finish()
     detected_blas(myenv)
     myenv.Append(LIBS=[blas_lib])
     return myenv
 else:
   print "Autodetecting BLAS support. See config.log for details!"
     
   ########################
   # check for acml
   myenv = env.Copy()
   # new configure object
   conf = Configure(myenv)
   # additional libs needed for acml
   myenv.Append(LIBS=['m', 'g2c'])
   if conf.CheckLib('acml', symname):
     detected_blas(myenv)
     found_blas = 1
   myenv = conf.Finish()
   if(found_blas == 1):
     return myenv
  
   ########################
   # check for goto
   myenv = env.Copy()
   # new configure object
   conf = Configure(myenv)
   # additional libs needed for goto
   myenv.Append(LIBS=['pthread'])
   # extra linker flags for libgoto
   # myenv.Append(_LIBFLAGS=['libs/numeric/mtl/build/xerbla.c'])
   myenv.Append(_LIBFLAGS=['libs/numeric/mtl/build/xerbla.o']) # not portable !!! 
   # myenv.Library('build/xerbla', 'build/xerbla.c')
   # myenv.Append(LIBS=['xerbla'])
   # myenv.Append(LIBPATH=['build'])

   if conf.CheckLib('goto', symname) or conf.CheckLib('goto_opteron-64', symname) or \
      conf.CheckLib('goto_coppermine-32', symname) or conf.CheckLib('goto_opteron-32', symname) or \
      conf.CheckLib('goto_itanium2-64', symname) or conf.CheckLib('goto_katmai-32', symname) or \
      conf.CheckLib('goto_northwood-32', symname) or conf.CheckLib('goto_prescott-32', symname) or \
      conf.CheckLib('goto_prescott-64', symname):
     detected_blas(myenv)
     found_blas = 1
   myenv = conf.Finish()
   if(found_blas == 1):
     return myenv
  
   ########################
   # check for ATLAS
   myenv = env.Copy()
   # new configure object
   conf = Configure(myenv)
   # additional libs needed for goto
   myenv.Append(LIBS=['f77blas', 'g2c'])
   if conf.CheckLib('atlas', symname):
     detected_blas(myenv)
     found_blas = 1
   myenv = conf.Finish()
   if(found_blas == 1):
     return myenv

 return env
 

 
###################################
# Add UMFPACK (extremely simplistic)
###################################

def detected_umfpack(env): # umfpack library found
   env.Append(CPPDEFINES = 'MTL_HAS_UMFPACK')

def check_for_umfpack(env):
 ufconfig_path = ARGUMENTS.get('ufconfig_path', '')
 if ufconfig_path:
   env.Prepend(CPPPATH = [ ufconfig_path ])
 else:
   print "UMFPACK only works with the UFconfig library."
   print "Provide ufconfig_path."
   print "Building continued without UMFPACK."
   return env

 amd_path = ARGUMENTS.get('amd_path', '')
 if amd_path:
   # print 'adding ' + amd_path + '/Lib to LIBPATH.'
   env.Append(LIBPATH = [ amd_path + '/Lib' ])
   #print 'adding ' + amd_path + '/Include to CPPPATH.'
   env.Prepend(CPPPATH = [ amd_path + '/Include' ])
   env.Prepend(LIBS=['amd'])
   # print "nachher: ", env['CPPPATH']
 else:
   print "UMFPACK only works with the AMD library."
   print "Provide amd_path."
   print "Building continued without UMFPACK."
   return env

 # get command line umfpack path
 umfpack_path = ARGUMENTS.get('umfpack_path', '')
 if umfpack_path:
   # print 'adding ' + umfpack_path + '/Lib to LIBPATH.'
   env.Prepend(LIBPATH = [ umfpack_path + '/Lib' ])
   env.Prepend(LIBS=['umfpack'])
   # print 'adding ' + umfpack_path + '/Include to CPPPATH.'
   env.Append(CPPPATH = [ umfpack_path + '/Include' ])
   detected_umfpack(env)
   return env
 



######################
# Start setting up env
######################


SourceSignatures('timestamp')

env = Environment()
# env.Decider('timestamp-newer')

conf = Configure(env)
check_no_long_double(conf)
env = conf.Finish()

# Search only for BLAS when explicitly asked for
with_blas = ARGUMENTS.get('with_blas', 0)
if int(with_blas):
    env = check_for_blas(env)

# Search only for UMFPACK when explicitly asked for
with_umfpack = ARGUMENTS.get('with_umfpack', 0)
if int(with_umfpack):
    env = check_for_umfpack(env)




######################
# Include paths
######################


my_cpppath = []
if os.environ.has_key('MTL_BOOST_ROOT') :
    mtlp = os.environ['MTL_BOOST_ROOT']
    my_cpppath.append(mtlp)
if os.environ.has_key('BOOST_ROOT') :
    my_cpppath.append(os.environ['BOOST_ROOT'])

# In branches, the current directory is prepended
pwd = os.getcwd()
if not os.environ.has_key('MTL_BOOST_ROOT') \
or not hasattr(os.path, 'samefile') or not os.path.samefile(pwd, mtlp):
     my_cpppath = [pwd] + my_cpppath

# Add directory of extern libraries if exist
# Not needed currently
#if os.path.exists('extern'):
#    my_cpppath.append(os.path.abspath('extern'))

env.Append(CPPPATH = my_cpppath)

# Hack to select compiler explicitly
# env['CXX'] = 'g++-4.1'


######################
# Opt. and debug flags
######################


opts = Options()
opts.Add('opt', 'Set to 1 for normal optimization, 2 for high optimization', 0)
opts.Add('debug', 'Set to 1 for debugging', 0)
opts.Add('check', 'Set to 1 if test programs shall be ran', 0)
opts.Add('with_concepts', 'Set to 1 if test programs shall be compiled with conceptg++ (must be in ~/bin)', 0)
opts.Add('full_warnings', 'Set to 1 to turn on all possible warnings and transform into errors (currently only gcc)', 0)


opts.Add('add_ccflag', 'Add extra CC compiler flags that aren\'t automatically set', '')
opts.Add('add_cxxflag', 'Add extra CXX compiler flags that aren\'t automatically set', '')
opts.Add('add_optflag', 'Add extra optimization flags that aren\'t automatically set', '')

############
# BLAS flags
############

opts.Add('with_blas', 'Link with BLAS (have to given each time :-( until we can save configurations)', 0)
opts.Add('blas_path', 'Add path of blas library', '')
opts.Add('blas_lib', 'Add name of blas library', '')
opts.Add('blas_ldflags', 'Add LDFLAGS for blas library', '')

############
# UMFPACK flags
############
opts.Add('with_umfpack', 'Use UMFPACK. Needs umfpack_path and amd_path.', 0)
opts.Add('umfpack_path', 'Add path of UMFPACK library', '')
opts.Add('amd_path', 'Add path of AMD library', '')
opts.Add('ufconfig_path', 'Add path of AMD library', '')


Help(opts.GenerateHelpText(env))

#env.Append(CCFLAGS = ARGUMENTS.get('debug', ''))
# Hack, for test only
# env.Replace(CXX = '/san/intel/cc/9.0/bin/icc')

# To test with conceptg++
# requires
if int(ARGUMENTS.get('with_concepts', 0)):
    env['CXX'] = '~/bin/conceptg++'
    env.Append(CCFLAGS = macro_flag(env, 'CONCEPTS_WITHOUT_OVERLOADED_REQUIREMENTS'))


# Full warnings and treat as errors
if int(ARGUMENTS.get('full_warnings', 0)):
    if env['CXX'] == 'g++':
        env.Append(CXXFLAGS = '-Werror -Wall -Wreorder -pedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings')
full_warnings = ARGUMENTS.get('full_warnings', 0)


# add user-defined CC flags 
add_ccflag = ARGUMENTS.get('add_ccflag', '')
if add_ccflag:
     env.Append(CCFLAGS = add_ccflag)

# add user-defined CXX flags 
add_cxxflag = ARGUMENTS.get('add_cxxflag', '')
if add_cxxflag:
     env.Append(CXXFLAGS = add_cxxflag)

add_platform_flags(env)

#create alternative environments
debug_env = env.Copy()
opt_env = env.Copy()
high_opt_env = env.Copy()

# add debugging flags to appropriate environment
debug = ARGUMENTS.get('debug', 0)
if int(debug):
    add_debug(env)
add_debug(debug_env)


# add optimization flags to appropriate environment
opt = ARGUMENTS.get('opt', 0)

# Add optimization flags and then copy C flags into C++ flags
add_opt(debug_env, 0) 
add_opt(env, int(opt))
add_opt(opt_env, 1)
add_opt(high_opt_env, 2)

# add user-defined optimization flags 
add_optflag = ARGUMENTS.get('add_optflag', '')
if add_optflag:
     add_user_opt(env, add_optflag)
     add_user_opt(opt_env, add_optflag)
     add_user_opt(high_opt_env, add_optflag)

# whether test programs should be ran
check = ARGUMENTS.get('check', 0)

#####################
# Sub-directories
#####################

Export('env debug_env opt_env high_opt_env check full_warnings')

SConscript(['libs/numeric/mtl/build/SConscript', 
            'libs/numeric/mtl/test/SConscript', 
            'libs/numeric/mtl/test_with_optimization/SConscript', 
            'libs/numeric/mtl/examples/SConscript', 
            'libs/numeric/mtl/experimental/SConscript',
            'libs/numeric/mtl/timing/SConscript', 
            'libs/numeric/linear_algebra/test/SConscript', 
            'libs/numeric/itl/test/SConscript'])

