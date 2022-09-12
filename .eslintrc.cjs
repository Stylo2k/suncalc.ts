module.exports = {
    extends: ['eslint:recommended', 'plugin:@typescript-eslint/recommended'],
    parser: '@typescript-eslint/parser',
    plugins: ['@typescript-eslint'],
    "parserOptions": {
        "sourceType": "module"
    },
    root: true,
    "rules": {
      "indent": 0,
      "array-bracket-spacing": 0,
      "strict": 0,
      "brace-style": 0
    },
    "env": {
      "amd": true
    }
};