/*=============================================================================
    Copyright (c) 2001-2007 Joel de Guzman
    Copyright (c) 2001-2009 Hartmut Kaiser

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
///////////////////////////////////////////////////////////////////////////////
//
//  A calculator example demonstrating the grammar and semantic actions
//  using phoenix to do the actual expression evaluation. The parser is
//  essentially an "interpreter" that evaluates expressions on the fly.
//
//  Additionally this examples shows how to build and use a lexer based on 
//  Ben Hansons Lexertl (http://www.benhanson.net/lexertl.html). This way the
//  parser matches the grammar against the tokens generated by the lexer 
//  component and not against the input character stream.
//
//  Even if the benefits of using a lexer for this small calculator grammar may 
//  not outweight the corresponding overhead, we provide this example because 
//  it allows to concentrate on the essentials without having to understand
//  the semantics first.
//
//  [ JDG June 29, 2002 ]   spirit1
//  [ JDG March 5, 2007 ]   spirit2
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/lex_lexer_lexertl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <iostream>
#include <string>

using namespace boost::spirit;
using namespace boost::spirit::qi;
using namespace boost::spirit::lex;
using namespace boost::spirit::ascii;
using namespace boost::spirit::arg_names;

///////////////////////////////////////////////////////////////////////////////
//  Our token definition
//  This class is used to define all the tokens to be recognized by the lexer. 
///////////////////////////////////////////////////////////////////////////////
template <typename Lexer>
struct calculator_tokens : lexer_def<Lexer>
{
    template <typename Self>
    void def (Self& self)
    {
        // unsigned integer token definition
        ui = "[1-9][0-9]*";
        
        // whitespace token definitions
        ws = "[ \\t\\f\\v]+";
        c_comment = "\\/\\*[^*]*\\*+([^/*][^*]*\\*+)*\\/";
        
        // build token set
        skipper = ws | c_comment;                   // += is allowed as well
        
        // associate the tokens and the token set with the lexer
        // default lexer state
        self = token_def<>('+') | '-' | '*' | '/' | '(' | ')'; 
        self += ui;                                 // still default state
        
        // The token_set 'skipper' get's assigned to a separate lexer state
        // which allows to use it separately from the main tokenization
        // (it is used as the skipper parser below)
        self("SKIPPER") = skipper;                  // lexer state "SKIPPER"
    }

    // This are the tokens to be recognized by the lexer.
    token_def<unsigned int> ui;   // matched tokens will have a unsigned int 
    token_def<> ws, c_comment;    // attribute will not be used

    // This is the only token set explicitly defined by this lexer because it
    // needs to be accessible from the outside (used as skip parser below).
    typename Lexer::token_set skipper;
};

///////////////////////////////////////////////////////////////////////////////
//  Our calculator grammar
//
//  The difference to the original example (calc3.cpp) is that we are 
//  specifying a second template parameter referring to the lexer. Further, we 
//  use a defined tokenset from above as the skip parser.
///////////////////////////////////////////////////////////////////////////////
template <typename Iterator, typename Lexer>
struct calculator : grammar<Iterator, int(), typename Lexer::token_set>
{
    template <typename TokenDef>
    calculator(TokenDef const& tok) 
      : calculator::base_type(expression)
    {
        // grammar
        expression =
            term                            [_val = _1]
            >> *(   ('+' >> term            [_val += _1])
                |   ('-' >> term            [_val -= _1])
                )
            ;

        term =
            factor                          [_val = _1]
            >> *(   ('*' >> factor          [_val *= _1])
                |   ('/' >> factor          [_val /= _1])
                )
            ;

        factor =
            tok.ui                          [_val = _1]
            |   '(' >> expression           [_val = _1] >> ')'
            |   ('-' >> factor              [_val = -_1])
            |   ('+' >> factor              [_val = _1])
            ;
    }

    rule<Iterator, int(), typename Lexer::token_set> expression, term, factor;
};

///////////////////////////////////////////////////////////////////////////////
//  Main program
///////////////////////////////////////////////////////////////////////////////
int
main()
{
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Expression parser...\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Type an expression...or [q or Q] to quit\n\n";

    // iterator type used to expose the underlying input stream
    typedef std::string::const_iterator base_iterator_type;
    
    // This is the lexer token type to use. The second template parameter lists 
    // all attribute types used for token_def's during token definition (see 
    // calculator_tokens<> above). Here we use the predefined lexertl token 
    // type, but any compatible token type may be used.
    typedef lexertl_token<
        base_iterator_type, boost::mpl::vector<unsigned int> 
    > token_type;
    
    // This is the lexer type to use to tokenize the input.
    // Here we use the lexertl based lexer engine.
    typedef lexertl_lexer<token_type> lexer_type;
    
    // This is the token definition type (derived from the given lexer type).
    typedef calculator_tokens<lexer_type> calculator_tokens;
    
    // this is the iterator type exposed by the lexer 
    typedef lexer<calculator_tokens>::iterator_type iterator_type;

    // this is the type of the grammar to parse
    typedef calculator<iterator_type, lexer_type> calculator;

    // now we use the types defined above to create the lexer and grammar
    // object instances needed to invoke the parsing process
    calculator_tokens tokens;                       // Our token definition
    calculator calc(tokens);                        // Our grammar definition

    lexer<calculator_tokens> lex(tokens);           // Our lexer

    // get input line by line and feed the parser to evaluate the expressions
    // read in from the input
    std::string str;
    int result;
    while (std::getline(std::cin, str))
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;

        // At this point we generate the iterator pair used to expose the
        // tokenized input stream.
        iterator_type iter = lex.begin(str.begin(), str.end());
        iterator_type end = lex.end();
        
        // Parsing is done based on the the token stream, not the character 
        // stream read from the input.
        // Note, how we use the token_set defined above as the skip parser.
        bool r = phrase_parse(iter, end, calc, result, tokens.skipper);

        if (r && iter == end)
        {
            std::cout << "-------------------------\n";
            std::cout << "Parsing succeeded\n";
            std::cout << "result = " << result << std::endl;
            std::cout << "-------------------------\n";
        }
        else
        {
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << "-------------------------\n";
        }
    }

    std::cout << "Bye... :-) \n\n";
    return 0;
}


