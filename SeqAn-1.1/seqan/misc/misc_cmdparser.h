#ifndef SEQAN_MISC_CMDPARSER
#define SEQAN_MISC_CMDPARSER

#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//  TODO:
//      * support multiple option values
//      * support some more formating options
//      * support options concatenated with their -fVALUE == -f VALUE
//      * store/return error code (invalid argument, invalid option, etc.)
//      * support named arguments (e.g. <ARG1> -> <INPUT FILE>)
//////////////////////////////////////////////////////////////////////////////
template<typename TChar>
inline bool
_isDigit(TChar const c)
{
    return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
            (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}

template<typename TString>
inline bool
_isDouble(TString const s)
{
    bool _dot = true;
    unsigned l = length(s);
    unsigned i = 0;

    // skip leading sign
    if(value (s,i) == '-') ++i;
    while(i < l){
        if(!_isDigit(value(s,i))){
            if(value(s,i) == '.' && _dot){
                _dot = false;
            }else return false;
        }
        ++i;
    }
    return true;
}

template<typename TString>
inline bool
_isInt(TString const s)
{
    unsigned l = length(s);
    unsigned i = 0;
    // skip leading sign
    if (value(s,i) == '-') ++i;
    while(i < l){
        if(!_isDigit(value(s,i))) return false;
        ++i;
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////

struct OptionType{
    enum {
        Boolean = 1,
        String = 2,
        Int = 4,
        Double = 8,
        Mandatory = 16,
        Debug = 32
    };
};

//////////////////////////////////////////////////////////////////////////////

/**
.Class.CommandLineOption:
..cat:Miscellaneous
..summary:Stores information for a specific Commandline Option
..signature:CommandLineOption
*/
class CommandLineOption{
public:
    CharString _longName;
    char         _shortName;

    CharString _helpText;
    int          _optionType;

    CommandLineOption() {}

    CommandLineOption(char _short,CharString _long,CharString _help,int _type)
        : _longName(_long),_shortName(_short),_helpText(_help),_optionType(_type)
        {
        }
/**.Memfunc.CommandLineOption#CommandLineOption:
..class:Class.CommandLineOption
..summary:Constructor
..signature:CommandLineOption ()
..signature:CommandLineOption (shortName,longName,helpText,type)
..param.shortName:A $char$ containing the short option identifier (e.g. $'h'$ for the $-h/--help$ option).
...remarks:Note that the leading "-" is not passed.
..param.longName:A @Shortcut.CharString@ containing the long option identifier (e.g. $"help"$ for the $-h/--help$ option).
...type:Shortcut.CharString
...remarks:Note that the leading "--" is not passed.
..param.helpText:A @Shortcut.CharString@ containing the help text associated with this option.
...type:Shortcut.CharString
*/
};

//////////////////////////////////////////////////////////////////////////////

/**
.Function.longName:
..summary:Returns the long option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:longName(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A @Shortcut.CharString@ holding the long name of the CommandLine Option (e.g. $help$ in case of $-h/--help$)
..remarks:The result type is @Shortcut.CharString@.
*/
inline CharString &
longName(CommandLineOption & me){
    return me._longName;
}

inline const CharString &
longName(CommandLineOption const & me){
    return me._longName;
}

/**
.Function.setLongName:
..summary:Sets the long option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setLongName(option,newName)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newName:A @Shortcut.CharString@ containing the new long name of the option.
...type:Shortcut.CharString
*/
inline void
setLongName(CommandLineOption & me, CharString const & newName){
    me._longName = newName;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.shortName:
..summary:Returns the short option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:shortName(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A $char$ holding the short name of the CommandLine Option (e.g. $h$ in case of $-h/--help$)
..remarks:The result type is $char$.
*/
inline char &
shortName(CommandLineOption & me){
    return me._shortName;
}

inline const char &
shortName(CommandLineOption const & me){
    return me._shortName;
}

/**
.Function.setShortName:
..summary:Sets the short option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setShortName(option,newName)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newName:A $char$ containing the new short name of the option.
*/
inline void
setShortName(CommandLineOption & me, char const & newName){
    me._shortName = newName;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.helpText:
..summary:Returns the help text associated with the @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:helpText(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A @Shortcut.CharString@ holding the help text of the CommandLine Option
..remarks:The result type is @Shortcut.CharString@.
*/
inline CharString &
helpText(CommandLineOption & me){
    return me._helpText;
}

inline const CharString &
helpText(CommandLineOption const & me){
    return me._helpText;
}

/**
.Function.setHelpText:
..summary:Sets the help text associated with the @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setHelpText(option,newHelpText)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newHelpText:A @Shortcut.CharString@ containing the new help text.
...type:Shortcut.CharString
*/
inline void
setHelpText(CommandLineOption & me, CharString const & newHelp){
    me._helpText = newHelp;
}

//////////////////////////////////////////////////////////////////////////////

inline const bool
isStringOption(CommandLineOption const & me){
    return ((me._optionType & OptionType::String) != 0);
}

inline const bool
isBooleanOption(CommandLineOption const & me){
    return ((me._optionType & OptionType::Boolean) != 0);
}

inline const bool
isDoubleOption(CommandLineOption const & me){
    return ((me._optionType & OptionType::Double) != 0);
}

inline const bool
isIntOption(CommandLineOption const & me){
    return ((me._optionType & OptionType::Int) != 0);
}

inline const bool
isDebugOption(CommandLineOption const & me){
    return ((me._optionType & OptionType::Debug) != 0);
}

inline const bool
isOptionMandatory(CommandLineOption const & me){
    return ((me._optionType & OptionType::Mandatory) != 0);
}

inline void
setOptionType(CommandLineOption & me,const int _newOpt){
    me._optionType = _newOpt;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream>
inline void
_writeOptName(TStream & target, CommandLineOption const & me)
{
    _streamWrite(target, ( shortName(me) == ' ' ? "" : "-"  ));
    if(shortName(me) != ' ') _streamPut(target,shortName(me));
    _streamWrite(target, ( shortName(me) == ' ' ||  longName(me) == "" ? "" : ", "  ));
    if(longName(me) != "")
    {
        _streamWrite(target, "--");
        _streamWrite(target, longName(me));
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream>
inline void
write(TStream & target, CommandLineOption const & me){
    _streamPut(target,'\t');
    _writeOptName(target, me);
    _streamPut(target,'\t');
    _streamPut(target,'\t');
    _streamWrite(target,me._helpText);
}

template <typename TStream>
inline TStream &
operator << (TStream & target, 
             CommandLineOption const & source)
{
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Class.CommandLineParser:
..cat:Miscellaneous
..summary:Stores multiple @Class.CommandLineOption@ objects and parses the command line arguments for these options
..signature:CommandLineParser
*/
class CommandLineParser{
public:
    typedef String<CommandLineOption>           TOptionMap;
    typedef Size<TOptionMap>::Type              TSize;
    
    typedef ::std::map<CharString, TSize>       TStringMap;
    typedef ::std::map<char, TSize >            TCharMap;
    typedef String < CharString >               TValueMap;

    String<CharString >  _cmdLine;
    TStringMap           _longNameMap;
    TCharMap             _shortNameMap;
    TValueMap            _valueMap;
    TOptionMap           _optionMap;
    
    unsigned             _required_arguments;
    String<CharString >  _arguments;
    CharString           _appName;


    unsigned line_width;
    unsigned padding_left;

/**.Memfunc.CommandLineParser#CommandLineParser:
..class:Class.CommandLineParser
..summary:Constructor
..signature:CommandLineParser ()
..signature:CommandLineParser (applicationName)
..param.applicationName:A @Shortcut.CharString@ containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

    CommandLineParser()
        : _required_arguments(0)
        {
            CommandLineOption opt('h',"help","displays this help message",OptionType::Boolean);
            appendValue(_optionMap,opt);
            insert(_shortNameMap,'h',0);
            insert(_longNameMap,"help",0);
            //insert(_long2ShortMap,longName(opt),shortName(opt));

            line_width   = 32;
            padding_left = 8;

            _appName = "";
        }

    CommandLineParser(CharString appName)
        : _required_arguments(0),_appName(appName)
        {
            CommandLineOption opt('h',"help","displays this help message",OptionType::Boolean);
            appendValue(_optionMap,opt);
            insert(_shortNameMap,'h',0);
            insert(_longNameMap,"help",0);

            line_width   = 32;
            padding_left = 8;
        }
};

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addOption:
..summary:adds an instance of @Class.CommandLineOption@ to the @Class.CommandLineParser@
..cat:Miscellaneous
..signature:addOption(parser,option)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The new @Class.CommandLineOption@ object that should be added.
...type:Class.CommandLineOption
*/
inline void
addOption(CommandLineParser & me,CommandLineOption const & opt){
    appendValue(me._optionMap,opt);
    if(shortName(opt) != ' ')  insert(me._shortNameMap,shortName(opt),length(me._optionMap) - 1);
    if(longName(opt) != "")  insert(me._longNameMap,longName(opt),length(me._optionMap) - 1);

    if(length(me._optionMap) > length(me._valueMap)) resize(me._valueMap , length(me._optionMap), Generous());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hasOption:
..summary:Returns true if the there is an option registered in the parser, that has the passed optionIdentifier
..cat:Miscellaneous
..signature:hasOption(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A $char$ or a @Shortcut.CharString@ that identifies the option.
*/
bool 
hasOption(CommandLineParser & me, CharString const & _long)
{
    return hasKey(me._longNameMap,_long);
}

bool 
hasOption(CommandLineParser const & me, CharString const & _long)
{
    return hasKey(me._longNameMap,_long);
}


bool 
hasOption(CommandLineParser & me, const char _short)
{
    return hasKey(me._shortNameMap,_short);
}

bool 
hasOption(CommandLineParser const & me, const char _short)
{
    return hasKey(me._shortNameMap,_short);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendCmdLine:
..summary:adds a line of text to the help output of the @Class.CommandLineParser@
..cat:Miscellaneous
..signature:appendCmdLine(parser,text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A text line that will be added to the help output.
...type:Shortcut.CharString
*/
inline void
appendCmdLine(CommandLineParser & me,CharString const & new_cmdLine){
    appendValue(me._cmdLine,new_cmdLine);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.requiredArguments:
..summary:using this option you can define how many non parameterized options are required by your program.
..cat:Miscellaneous
..signature:requireRemainder(parser,count)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.count:A $unsigned int$ defining the amount of non-parameterized options requried by your program.
*/
inline void
requiredArguments(CommandLineParser & me,unsigned count){
    me._required_arguments = count;
}	


//////////////////////////////////////////////////////////////////////////////

template <typename TStream>
inline void
_usage(CommandLineParser & me, TStream & target)
{
    _streamWrite(target, "Usage: ");
    _streamWrite(target, me._appName);
    _streamWrite(target, " [OPTION]... ");
    for(unsigned r = 0; r < me._required_arguments; ++r)
    {
        _streamWrite(target, "<ARG");
        _streamPutInt(target, r+1);
        _streamWrite(target,"> ");
    }
    _streamPut(target,'\n');
}

/**
.Function.shortHelp:
..summary:Prints a short help message for your parser to the stream
..cat:Miscellaneous
..signature:shortHelp(parser,stream)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
*/
template <typename TStream>
inline void
shortHelp(CommandLineParser & me, TStream & target){
    _usage(me,target);
    _streamWrite(target, "Try '");
    _streamWrite(target, me._appName);
    _streamWrite(target, " --help' for more information.\n");
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.help:
..summary:Prints the complete help message for your parser to the stream.
..cat:Miscellaneous
..signature:help(parser[,stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
*/
template <typename TStream>
inline void
help(CommandLineParser & me, TStream & target){
    //_streamWrite(target,"Usage: ");
    for(unsigned i = 0; i < length(me._cmdLine);++i){ _streamWrite(target,me._cmdLine[i]);_streamPut(target,'\n');}

    _streamPut(target,'\n');
    _usage(me,target);
    _streamPut(target,'\n');

    for(unsigned o = 0;o < length(me._optionMap);++o)
    {
        const CommandLineOption opt = value(me._optionMap,o);       
        if(isDebugOption(opt)) continue;    // do not print debug options .. these are not for the user
        unsigned s = 0;
       
        while(s < me.padding_left){
            _streamPut(target,' ');
            ++s;
        }

        _streamPut(target, ( shortName(opt) == ' ' ? ' ' : '-'  ));
        _streamPut(target,shortName(opt));
        
        s += 2;
        
        _streamPut(target, ( shortName(opt) == ' ' ||  longName(opt) == "" ? ' ' : ','  ));++s;
        
        if(longName(opt) != "")
        {
            s += 3; // 4 signs => '-x, '
            _streamWrite(target, " --");
            _streamWrite(target, longName(opt));
            s += length(longName(opt));
        }

        if(s < me.line_width){
            while(s < me.line_width)
            {
                _streamPut(target,' ');
                ++s;
            }
            _streamWrite(target,helpText(opt));
        }
        else
        {
            _streamPut(target,'\n');
            s = 0;
            while(s < me.line_width){
                _streamPut(target,' ');
                ++s;
            }
            _streamWrite(target,helpText(opt));
        }
        _streamPut(target,'\n');
    }
}

inline void
help(CommandLineParser & me)
{
    help(me,::std::cerr);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.isSet:
..summary:Returns true if the option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A $char$ or a @Shortcut.CharString@ that identifies the option.
*/
inline bool
isSet(CommandLineParser & me,char const & shortName)
{
    if(!hasKey(me._shortNameMap,shortName)) return false; // this option does not exist
    else
    {
        // if value != "" -> value was set
        if(value(me._valueMap,cargo(me._shortNameMap,shortName)) != "") return true;
        else return false;
    }
}

inline bool
isSet(CommandLineParser & me,CharString const & longName)
{
    if(!hasKey(me._longNameMap,longName)) return false; // this option does not exist
    else
    {
        // if value != "" -> value was set
        if(value(me._valueMap,cargo(me._longNameMap,longName)) != "") return true;
        else return false;
    }
}

//////////////////////////////////////////////////////////////////////////////

inline bool
_allMandatorySet(CommandLineParser & me)
{
    for(unsigned o = 0;o < length(me._optionMap);++o)
        if(value(me._valueMap,o) == "" && isOptionMandatory(value(me._optionMap,o))) return false;
    return true;
}

//////////////////////////////////////////////////////////////////////////////

inline CharString
_parseAppName(CharString const & candidate)
{
    int l = length(candidate);
    int i;

    for(i = l - 1; l >= 0;--i)
        if(value(candidate,i) == '\\' || value(candidate,i) == '/') 
            break;
    ++i;
    CharString ret = "";
    for(int j = i;j < l;++j) append(ret,value(candidate,j));
    return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TErrorStream>
bool _assignOptionValue(CommandLineParser & me, unsigned option_index, CharString const & _val, TErrorStream & estream)
{
    // get the option object
    CommandLineOption opt = value(me._optionMap,option_index);
    if(isDoubleOption(opt)){
        if(!_isDouble(_val))
        {
            _streamWrite(estream,me._appName);
            _streamWrite(estream,": ");
            _streamWrite(estream, "\"");
            _streamWrite(estream, _val);
            _streamWrite(estream, "\" is not a valid double value for '");
            _writeOptName(estream, opt);
            _streamWrite(estream, "'\n");
            return false;
        }
    }else if(isIntOption(opt)){
        if(!_isInt(_val))
        {
            _streamWrite(estream,me._appName);
            _streamWrite(estream,": ");
            _streamWrite(estream, "\"");
            _streamWrite(estream, _val);
            _streamWrite(estream, "\" is not a valid integer value for '");
            _writeOptName(estream, opt);
            _streamWrite(estream, "'\n");
            return false;
        }
    }    
    value(me._valueMap,option_index) = _val;
    return true;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.parse:
..summary:Returns true if the option was set on the parsed command line.
..cat:Miscellaneous
..signature:parse(parser,argc,argv[,errorStream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.argc:Count of the objects on the command line.
..param.argv:Array of the different command line arguments ($const char *argv[]$). 
..param.errorStream:A stream where error messages are send too.
*/
template<typename TErrorStream>
bool
parse(CommandLineParser & me,int argc, const char *argv[], TErrorStream & estream)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;
    // if the appName wasn't set .. parse from command line
    if(me._appName == "") me._appName = _parseAppName(argv[0]);

    for(int i = 1; i < argc; ++i) 
    {
        if (argv[i][0] == '-')  // this is possibly an option value
        {
            CharString inParam = argv[i];
            unsigned len = length (inParam);
            char shortOpt = '-'; // 
            CharString longOpt;
            
            if(len == 1)
            {
                _streamWrite(estream,me._appName);
                _streamWrite(estream,": invalid option '-'\n");
                return false;
            }
            else if(len == 2)
            {
                if(value(inParam,1) != '-')
                {
                    shortOpt = value(inParam,1);
                    if(hasOption(me,shortOpt)) // this is a short option
                    {
                        TOptionPosition option_index = cargo(me._shortNameMap,shortOpt);
                        CommandLineOption opt = value(me._optionMap,option_index);

                        if(!isBooleanOption(opt) && (i + 1) == argc) // no value available
                        {
                            _streamWrite(estream,me._appName);
                            _streamWrite(estream,": ");
                            _streamWrite(estream, "'");
                            _writeOptName(estream, opt);
                            _streamWrite(estream, "' requires a value\n");
                            return false;
                        }
                        else if(isBooleanOption(opt))
                        {
                            value(me._valueMap,option_index) = "set";
                        }
                        else
                        {
                            ++i;
                            CharString _val = argv[i];
                            if(!_assignOptionValue(me,option_index,_val,estream)) return false;
                        }
                    }else{ // ERROR -> the parser does not recognize this option
                        _streamWrite(estream,me._appName);
                        _streamWrite(estream,": invalid option ");
                        _streamWrite(estream, "'-");
                        _streamPut(estream, shortOpt);
                        _streamWrite(estream, "'\n");
                        return false;
                    }
                }else{
                    _streamWrite(estream,me._appName);
                    _streamWrite(estream,": invalid option '--'\n");
                    return false;
                }
            }
            else if (value(inParam,0) == '-' && value(inParam,1) != '-') // maybe a combination of multiple bool opts
            {
                for(unsigned o = 1; o < len;++o)
                {
                    shortOpt = value(inParam,o);
                    if(!hasOption(me,shortOpt))
                    {
                        _streamWrite(estream,me._appName);
                        _streamWrite(estream,": invalid option ");
                        _streamWrite(estream, "'-");
                        _streamPut(estream, shortOpt);
                        _streamWrite(estream, "'\n");
                        return false;
                    }
                        TOptionPosition option_index = cargo(me._shortNameMap,shortOpt);
                    CommandLineOption opt = value(me._optionMap,option_index);
                    if(!isBooleanOption(opt)) // value options can not be part of multi options
                    {
                        _streamWrite(estream,me._appName);
                        _streamWrite(estream,": ");
                        _streamPut(estream, '\'');
                        _writeOptName(estream, opt);
                        _streamWrite(estream, "' requires a value and therefore can not be part of the argument list ");
                        _streamWrite(estream, inParam);
                        _streamPut(estream, '\n');
                        return false;                        
                    }else
                        value(me._valueMap,option_index) = "set";
                }
            }
            else if (value(inParam,0) == '-' && value(inParam,1) == '-') // this is a long option
            {
                unsigned t = 2;
                CharString _val;
                while(t < len && value(inParam,t) != '=')
                {
                    appendValue(longOpt,value(inParam,t));
                    ++t;
                }
                if(t < len) // this one is a --name=value option
                {
                    ++t;
                    while(t < len)
                    {
                        appendValue(_val,value(inParam,t));
                        ++t;
                    }
                }
                // we may be got already a value
                if(hasOption(me,longOpt))
                {
                    TOptionPosition option_index = cargo(me._longNameMap,longOpt);
                    CommandLineOption opt = value(me._optionMap,option_index);

                    if ( _val != "")
                    {
                        if(!_assignOptionValue(me,option_index,_val,estream)) return false;
                    }
                    else if((!isBooleanOption(opt)) && ((i + 1) == argc)) // no value available
                    {
                        _streamWrite(estream,me._appName);
                        _streamWrite(estream,": ");
                        _streamWrite(estream, "'");
                        _writeOptName(estream, opt);
                        _streamWrite(estream, "' requires a value\n");
                        return false;
                    }
                    else if(isBooleanOption(opt))
                    {
                        value(me._valueMap,option_index) = "set";
                    }
                    else
                    {
                        ++i;
                        _val = argv[i];
                        if(!_assignOptionValue(me,option_index,_val,estream)) return false;
                    }
                }
                else
                {
                    _streamWrite(estream,me._appName);
                    _streamWrite(estream,": invalid option");
                    _streamWrite(estream, "'--");
                    _streamWrite(estream, longOpt);
                    _streamWrite(estream, "'\n");
                    return false;
                }
            }            
        }
        else
        { // this seems to be a normal argument
            appendValue(me._arguments,argv[i] );
        }
    }
    if(isSet(me,'h'))
    {
        //help(me, estream);
        return true;
    }
    else if(!_allMandatorySet(me) || length(me._arguments) < me._required_arguments) 
        return false;
    else 
        return true;
}

//////////////////////////////////////////////////////////////////////////////

inline bool
parse(CommandLineParser & me,int argc, const char *argv[])
{
    return parse(me,argc,argv,::std::cerr);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.getOptionValue:
..summary:Fills the passed variable $value$ with the value set for the option on the command line.
..cat:Miscellaneous
..signature:getOptionValue(parser,optionIdentifier,value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A $char$ or a @Shortcut.CharString@ that identifies the option.
..param.value:The variable where the value is stored.
...remarks: The variable type ($int$, $double$, $bool$ or @Shortcut.CharString@) depends on the OptionTyp.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
*/
inline bool
getOptionValue(CommandLineParser & me,char const & shortName, bool & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,shortName)) return false;
    TOptionPosition option_index = cargo(me._shortNameMap,shortName);
    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isBooleanOption(opt)) return false;
    else{
        val = !(value(me._valueMap,option_index) == "");
        return true;
    }
}

inline bool
getOptionValue(CommandLineParser & me,CharString const & longName, bool & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,longName)) return false;
    TOptionPosition option_index = cargo(me._longNameMap,longName);
    
    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isBooleanOption(opt)) return false;
    else{
        val = !(value(me._valueMap,option_index) == "");
        return true;
    }
}

//////////////////////////////////////////////////////////////////////////////

inline bool
getOptionValue(CommandLineParser & me,char const & shortName, CharString & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,shortName)) return false;
    TOptionPosition option_index = cargo(me._shortNameMap,shortName);

    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isStringOption(opt)) return false;
    else{
        val = value(me._valueMap,option_index);
        return isSet(me,shortName);
    }
}

inline bool
getOptionValue(CommandLineParser & me,CharString const & longName, CharString & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,longName)) return false;
    TOptionPosition option_index = cargo(me._longNameMap,longName);
    
    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isStringOption(opt)) return false;
    else{
        val = value(me._valueMap,option_index);
        return isSet(me,longName);
    }
}

//////////////////////////////////////////////////////////////////////////////

inline bool
getOptionValue(CommandLineParser & me,char const & shortName, int & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,shortName)) return false;
    TOptionPosition option_index = cargo(me._shortNameMap,shortName);

    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isIntOption(opt)) return false;
    else{
        val = atoi(toCString(value(me._valueMap,option_index)));
        return isSet(me,shortName);
    }
}

inline bool
getOptionValue(CommandLineParser & me,CharString const & longName, int & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,longName)) return false;
    TOptionPosition option_index = cargo(me._longNameMap,longName);
    
    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isIntOption(opt)) return false;
    else{
        val = atoi(toCString(value(me._valueMap,option_index)));
        return isSet(me,longName);
    }
}

//////////////////////////////////////////////////////////////////////////////

inline bool
getOptionValue(CommandLineParser & me,char const & shortName, double & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,shortName)) return false;
    TOptionPosition option_index = cargo(me._shortNameMap,shortName);

    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isDoubleOption(opt)) return false;
    else{
        val = atof(toCString(value(me._valueMap,option_index)));
        return isSet(me,shortName);
    }
}

inline bool
getOptionValue(CommandLineParser & me,CharString const & longName, double & val){
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if(!hasOption(me,longName)) return false;
    TOptionPosition option_index = cargo(me._longNameMap,longName);
    
    CommandLineOption opt = value(me._optionMap,option_index);
    if(!isDoubleOption(opt)) return false;
    else{
        val = atof(toCString(value(me._valueMap,option_index)));
        return isSet(me,longName);
    }
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getArgumentValue:
..summary:Fills the passed variable $value$ with the argument set on the command line.
..cat:Miscellaneous
..signature:getArgumentValue(parser,position,value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.position:A zero based $int4 indicating which argument you want to get.
..param.value:The variable where the value is stored.
...type:Shortcut.CharString
..returns: $true$ if the requested argument exists, $false$ otherwise.
*/
inline bool
getArgumentValue(CommandLineParser & me, unsigned position, CharString & val){
    if(position < length(me._arguments)){
        val = me._arguments[position];
        return true;
    }else return false;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.argumentCount:
..summary:Returns the count of passed arguments.
..cat:Miscellaneous
..signature:argumentCount(parser)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
*/
inline Size<String<CharString> >::Type
argumentCount(CommandLineParser & me){
    return length(me._arguments);
}


} // end SEQAN_NAMESPACE_MAIN

#endif
