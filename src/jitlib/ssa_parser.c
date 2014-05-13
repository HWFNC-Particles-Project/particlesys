
#include "ssa_def.h"
#include <setjmp.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "hash_table.h"

void skip_space(char *str, int *i) {
    while(isspace(str[*i])) ++*i;
}

enum token_t {
    TOKEN_NONE,
    TOKEN_NUMBER,
    TOKEN_IDENTIFIER,
    TOKEN_LPAREN,
    TOKEN_RPAREN,
    TOKEN_LBRACKET,
    TOKEN_RBRACKET,
    TOKEN_COMMA,
    TOKEN_ASSIGN,
    TOKEN_PLUS,
    TOKEN_MINUS,
    TOKEN_MUL,
    TOKEN_DIV,
    TOKEN_QUESTION,
    TOKEN_COLON,
    TOKEN_EQ,
    TOKEN_LT,
    TOKEN_LE,
    TOKEN_NEQ,
    TOKEN_NLT,
    TOKEN_NLE,
    TOKEN_AND,
    TOKEN_ANDN,
    TOKEN_OR,
    TOKEN_XOR,
    TOKEN_SQRT,
    TOKEN_RSQRT,
    TOKEN_RCP,
    TOKEN_MIN,
    TOKEN_MAX,
};

typedef struct {
    char *str;
    int i, begin;
    float number;
    enum token_t current;
    int capacity, size;
    uint64_t *ssa_out;
    hash_table symbols;
    jmp_buf error_jmp;
} parser_state;

enum token_t next_token(parser_state *parser) {
    char *str = parser->str;
    int *i = &parser->i;
    skip_space(str, i);
    parser->begin = *i;
    if(str[*i] == '+') {
        ++*i; return TOKEN_PLUS;
    } else if(str[*i] == '-') {
        ++*i; return TOKEN_MINUS;
    } else if(str[*i] == '*') {
        ++*i; return TOKEN_MUL;
    } else if(str[*i] == '/') {
        ++*i; return TOKEN_DIV;
    } else if(str[*i] == '(') {
        ++*i; return TOKEN_LPAREN;
    } else if(str[*i] == ')') {
        ++*i; return TOKEN_RPAREN;
    } else if(str[*i] == '[') {
        ++*i; return TOKEN_LBRACKET;
    } else if(str[*i] == ']') {
        ++*i; return TOKEN_RBRACKET;
    } else if(str[*i] == ',') {
        ++*i; return TOKEN_COMMA;
    } else if(str[*i] == ':') {
        ++*i; return TOKEN_COLON;
    } else if(str[*i] == '?') {
        ++*i; return TOKEN_QUESTION;
    } else if(str[*i] == '|') {
        ++*i; return TOKEN_OR;
    } else if(str[*i] == '^') {
        ++*i; return TOKEN_XOR;
    } else if(str[*i] == '&') {
        ++*i;
        if(str[*i] == '!') {
            ++*i; return TOKEN_ANDN;
        } else {
            return TOKEN_AND;
        }
    } else if(str[*i] == '=') {
        ++*i;
        if(str[*i] == '=') {
            ++*i; return TOKEN_EQ;
        } else {
            return TOKEN_ASSIGN;
        }
    } else if(str[*i] == '<') {
        ++*i;
        if(str[*i] == '=') {
            ++*i; return TOKEN_LE;
        } else {
            return TOKEN_LT;
        }
    } else if(str[*i] == '>') {
        ++*i;
        if(str[*i] == '=') {
            ++*i; return TOKEN_NLT;
        } else {
            return TOKEN_NLE;
        }
    } else if(str[*i] == 's' && str[*i+1] == 'q' && str[*i+2] == 'r' && str[*i+3] == 't') {
        *i+=4; return TOKEN_SQRT;
    } else if(str[*i] == 'r' && str[*i+1] == 's' && str[*i+2] == 'q' && str[*i+3] == 'r' && str[*i+4] == 't') {
        *i+=5; return TOKEN_RSQRT;
    } else if(str[*i] == 'r' && str[*i+1] == 'c' && str[*i+2] == 'p') {
        *i+=3; return TOKEN_RCP;
    } else if(str[*i] == 'm' && str[*i+1] == 'i' && str[*i+2] == 'n') {
        *i+=3; return TOKEN_MIN;
    } else if(str[*i] == 'm' && str[*i+1] == 'a' && str[*i+2] == 'x') {
        *i+=3; return TOKEN_MAX;
    } else if(isalpha(str[*i])) {
        while(isalnum(str[*i])) ++*i;
        return TOKEN_IDENTIFIER;
    } else if(isdigit(str[*i]) || str[*i] == '.') {
        char *end;
        parser->number = strtof(str+*i, &end);
        *i = end-str;
        return TOKEN_NUMBER;
    } else if(str[*i] == '\0') {
        return TOKEN_NONE;
    } else {
        fprintf(stderr, "unexpected character: '%c'\n", str[*i]);
        longjmp(parser->error_jmp, 1);
        return TOKEN_NONE;
    }
}

parser_state parser_init(char *str) {
    parser_state parser;
    parser.str = str;
    parser.i = 0;
    parser.capacity = 4;
    parser.size = 0;
    parser.ssa_out = realloc(NULL, parser.capacity*sizeof(uint64_t));
    parser.current = next_token(&parser);
    parser.symbols = hash_table_create();
    return parser;
}

void parser_free(parser_state *parser) {
    hash_table_destroy(parser->symbols);
}

int parser_append_ssa(parser_state *parser, uint64_t ssa) {
    if(parser->size == parser->capacity) {
        parser->capacity = parser->capacity*3/2;
        parser->ssa_out = realloc(parser->ssa_out, parser->capacity*sizeof(uint64_t));
    }
    parser->ssa_out[parser->size] = ssa;
    return parser->size++;
}

void advance(parser_state *parser) {
    parser->current = next_token(parser);
}
int peek(parser_state *parser, enum token_t t) {
    return parser->current == t;
}
int accept(parser_state *parser, enum token_t t) {
    if(peek(parser, t)) {
        advance(parser);
        return 1;
    } else {
        return 0;
    }
}

int expect(parser_state *parser, enum token_t t) {
    if(accept(parser, t)) {
        return 1;
    } else {
        fprintf(stderr, "unexpected token: '%.*s'\n", parser->i-parser->begin, parser->str+parser->begin);
        longjmp(parser->error_jmp, 1);
        return 0;
    }
}

int expression(parser_state *parser);

int primary(parser_state *parser) {
    if(peek(parser, TOKEN_IDENTIFIER)) {
        int *value = hash_table_get_n(parser->symbols, parser->str+parser->begin, parser->i - parser->begin);
        if(value == NULL) {
            fprintf(stderr, "undefined identifier: %.*s\n", parser->i - parser->begin, parser->str+parser->begin);
            longjmp(parser->error_jmp, 1);
        }
        advance(parser);
        return *value;
    } else if(accept(parser, TOKEN_SQRT)) {
        expect(parser, TOKEN_LPAREN);
        int arg = expression(parser);
        int value = parser_append_ssa(parser, ssa_node(SSA_SQRT, arg, 0, 0));
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(accept(parser, TOKEN_RSQRT)) {
        expect(parser, TOKEN_LPAREN);
        int arg = expression(parser);
        int value = parser_append_ssa(parser, ssa_node(SSA_RSQRT, arg, 0, 0));
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(accept(parser, TOKEN_RCP)) {
        expect(parser, TOKEN_LPAREN);
        int arg = expression(parser);
        int value = parser_append_ssa(parser, ssa_node(SSA_RCP, arg, 0, 0));
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(accept(parser, TOKEN_MIN)) {
        expect(parser, TOKEN_LPAREN);
        int arg1 = expression(parser);
        expect(parser, TOKEN_COMMA);
        int arg2 = expression(parser);
        int value = parser_append_ssa(parser, ssa_node(SSA_MIN, arg1, arg2, 0));
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(accept(parser, TOKEN_MAX)) {
        expect(parser, TOKEN_LPAREN);
        int arg1 = expression(parser);
        expect(parser, TOKEN_COMMA);
        int arg2 = expression(parser);
        int value = parser_append_ssa(parser, ssa_node(SSA_MAX, arg1, arg2, 0));
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(peek(parser, TOKEN_NUMBER)) {
        float f = parser->number;
        int value = parser_append_ssa(parser, ssa_const(f));
        advance(parser);
        return value;
    } else if(accept(parser, TOKEN_LPAREN)) {
        int value = expression(parser);
        expect(parser, TOKEN_RPAREN);
        return value;
    } else if(accept(parser, TOKEN_LBRACKET)) {
        float f = parser->number;
        int value = parser_append_ssa(parser, ssa_node(SSA_LOAD, f, 0, 0));
        expect(parser, TOKEN_NUMBER);
        expect(parser, TOKEN_RBRACKET);
        return value;
    } else {
        fprintf(stderr, "expected primary expression\n");
        longjmp(parser->error_jmp, 1);
        return -1;
    }
}

int multiplicative(parser_state *parser) {
    int left = primary(parser);
    for(;;) {
        if(accept(parser, TOKEN_MUL)) {
            int right = primary(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_MUL, left, right, 0));
        } else if(accept(parser, TOKEN_DIV)) {
            int right = primary(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_DIV, left, right, 0));
        } else {
            return left;
        }
    }
}

int additive(parser_state *parser) {
    int left = multiplicative(parser);
    for(;;) {
        if(accept(parser, TOKEN_PLUS)) {
            int right = multiplicative(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_ADD, left, right, 0));
        } else if(accept(parser, TOKEN_MINUS)) {
            int right = multiplicative(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_SUB, left, right, 0));
        } else {
            return left;
        }
    }
}

int relational(parser_state *parser) {
    int left = additive(parser);
    for(;;) {
        if(accept(parser, TOKEN_LT)) {
            int right = additive(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPLT, left, right, 0));
        } else if(accept(parser, TOKEN_LE)) {
            int right = additive(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPLE, left, right, 0));
        } else if(accept(parser, TOKEN_NLT)) {
            int right = additive(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPNLT, left, right, 0));
        } else if(accept(parser, TOKEN_NLE)) {
            int right = additive(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPNLE, left, right, 0));
        } else {
            return left;
        }
    }
}

int equality(parser_state *parser) {
    int left = relational(parser);
    for(;;) {
        if(accept(parser, TOKEN_EQ)) {
            int right = relational(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPEQ, left, right, 0));
        } else if(accept(parser, TOKEN_NEQ)) {
            int right = relational(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_CMPNEQ, left, right, 0));
        } else {
            return left;
        }
    }
}

int and(parser_state *parser) {
    int left = equality(parser);
    for(;;) {
        if(accept(parser, TOKEN_AND)) {
            int right = equality(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_AND, left, right, 0));
        } else if(accept(parser, TOKEN_ANDN)) {
            int right = equality(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_ANDN, left, right, 0));
        } else {
            return left;
        }
    }
}

int xor(parser_state *parser) {
    int left = and(parser);
    for(;;) {
        if(accept(parser, TOKEN_XOR)) {
            int right = and(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_XOR, left, right, 0));
        } else {
            return left;
        }
    }
}

int or(parser_state *parser) {
    int left = xor(parser);
    for(;;) {
        if(accept(parser, TOKEN_OR)) {
            int right = xor(parser);
            left = parser_append_ssa(parser, ssa_node(SSA_OR, left, right, 0));
        } else {
            return left;
        }
    }
}

int ternary(parser_state *parser) {
    int left = or(parser);
    if(accept(parser, TOKEN_QUESTION)) {
        int middle = ternary(parser);
        expect(parser, TOKEN_COLON);
        int right = ternary(parser);
        return parser_append_ssa(parser, ssa_node(SSA_BLEND, left, middle, right));
    } else {
        return left;
    }
}

int expression(parser_state *parser) {
    return ternary(parser);
}

void ssa_block_load(ssa_block *block, uint64_t *code, size_t size);

ssa_block ssa_parse(char *str) {
    parser_state parser = parser_init(str);

    if(setjmp(parser.error_jmp)) {
        free(parser.ssa_out);
        parser_free(&parser);
        ssa_block block;
        ssa_block_load(&block, NULL, 0);
        return block;
    }

    while(!peek(&parser,TOKEN_NONE)) {
        if(peek(&parser, TOKEN_IDENTIFIER)) {
            int *value = hash_table_insert_n(parser.symbols, parser.str + parser.begin, parser.i - parser.begin);
            advance(&parser);
            expect(&parser, TOKEN_ASSIGN);
            *value = expression(&parser);
        } else if(accept(&parser, TOKEN_LBRACKET)) {
            int target = strtof(parser.str+parser.begin, NULL);
            expect(&parser, TOKEN_NUMBER);
            expect(&parser, TOKEN_RBRACKET);
            expect(&parser, TOKEN_ASSIGN);
            int value = expression(&parser);
            parser_append_ssa(&parser, ssa_node(SSA_STORE, value, target, 0));
        } else {
            fprintf(stderr, "unexpected token: '%.*s'\n", parser.i-parser.begin, parser.str+parser.begin);
            longjmp(parser.error_jmp, 1);
        }
    }
    ssa_block block;
    block.size = parser.size;
    block.capacity = parser.capacity;
    block.buffer = parser.ssa_out;
    block.remap = malloc(parser.capacity*sizeof(uint32_t));
    for(uint32_t i = 0;i<block.size;++i) {
        block.remap[i] = i;
    }
    parser_free(&parser);

    return block;
}
